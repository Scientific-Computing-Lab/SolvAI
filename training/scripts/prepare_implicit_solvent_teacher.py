"""Build leakage-safe teachers from public GBn2 and learned-solvent campaigns.

The source contains molecule-level free energies for nearly all FreeSolv
molecules from a GBn2 alchemical calculation and from a neural implicit-solvent
potential trained on 280k explicit-solvent frames.  Benchmark identities are
removed globally before fitting the structure-only teacher.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict

import joblib
import numpy as np
import pandas as pd
from arrow_distill.data import FREESOLV_JSON, ROOT, canonicalize
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold

SOURCE_DIR = ROOT / "data/external/repos/mlimplicitsolvent-RAW"
TARGETS = (
    "gbn2_alchemical_dg",
    "neural_implicit_dg",
    "gbn2_exp_residual",
    "neural_implicit_exp_residual",
    "neural_minus_gbn2",
)


def normalize_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def estimator(seed: int, leaf: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=500,
        max_features=0.7,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )


def load_source(benchmark_keys: set[str]) -> pd.DataFrame:
    freesolv = json.loads(FREESOLV_JSON.read_text())
    names: dict[str, set[str]] = defaultdict(set)
    for public_id, record in freesolv.items():
        for name in (record.get("iupac"), record.get("nickname")):
            if name:
                names[normalize_name(name)].add(public_id)

    explicit = pd.read_csv(SOURCE_DIR / "explicit.csv")
    neural = pd.read_csv(SOURCE_DIR / "yank.csv").rename(
        columns={"Molecule": "iupac", " Implicit Free Energy of Solvation": "neural_implicit_dg"}
    )
    neural_by_name = neural.groupby(neural.iupac.map(normalize_name)).neural_implicit_dg.median()
    rows: list[dict[str, object]] = []
    for record in explicit.itertuples(index=False):
        candidates = names[normalize_name(record.iupac)]
        if not candidates:
            continue
        # The only ambiguous name is a stereoisomeric pair with one connectivity.
        public_id = min(candidates)
        source = freesolv[public_id]
        smiles, inchi_key, connectivity_key, _ = canonicalize(source["smiles"])
        if connectivity_key in benchmark_keys:
            continue
        neural_value = neural_by_name.get(normalize_name(record.iupac), np.nan)
        rows.append(
            {
                "public_id": ";".join(sorted(candidates)),
                "canonical_smiles": smiles,
                "inchi_key": inchi_key,
                "connectivity_key": connectivity_key,
                "gbn2_alchemical_dg": float(record.delta_F),
                "neural_implicit_dg": float(neural_value),
                "gbn2_exp_residual": float(record.exp_dG - record.delta_F),
                "neural_implicit_exp_residual": float(record.exp_dG - neural_value),
                "neural_minus_gbn2": float(neural_value - record.delta_F),
            }
        )
    frame = pd.DataFrame(rows)
    aggregations: dict[str, object] = {target: "median" for target in TARGETS}
    aggregations.update(
        {
            "public_id": lambda values: ";".join(sorted(set(values))),
            "canonical_smiles": "first",
            "inchi_key": "first",
        }
    )
    return frame.groupby("connectivity_key", as_index=False).agg(aggregations)


def main() -> None:
    processed = ROOT / "data/processed"
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    benchmark_features = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    expanded_path = processed / "expanded_public_hydration_nonbenchmark.parquet"
    expanded_features_path = processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    if expanded_path.exists() and expanded_features_path.exists():
        public = pd.read_parquet(expanded_path)
        public_features = pd.read_parquet(expanded_features_path)
    else:
        public = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
        public_features = pd.read_parquet(
            processed / "public_hydration_rdkit_morgan_features.parquet"
        )
    benchmark_keys = set(benchmark.inchi_connectivity_key)
    source = load_source(benchmark_keys)
    if set(source.connectivity_key) & benchmark_keys:
        raise AssertionError("ARROW benchmark connectivity leaked into implicit teacher")
    source.to_parquet(processed / "implicit_solvent_nonbenchmark.parquet", index=False)

    descriptor_columns = [
        column for column in benchmark_features if column.startswith(("rdkit__", "morgan2__"))
    ]
    x_benchmark = benchmark_features[descriptor_columns].to_numpy(dtype=np.float32)
    x_public = public_features[descriptor_columns].to_numpy(dtype=np.float32)
    public_position = {
        key: index for index, key in enumerate(public.inchi_connectivity_key.astype(str))
    }
    benchmark_output = benchmark[
        ["molecule_id", "canonical_smiles", "inchi_connectivity_key"]
    ].rename(columns={"inchi_connectivity_key": "connectivity_key"})
    public_output = public[["molecule_id", "canonical_smiles", "inchi_connectivity_key"]].rename(
        columns={"inchi_connectivity_key": "connectivity_key"}
    )
    benchmark_output["prediction_scope"] = "benchmark_external_teacher"
    scope_has_oof = np.zeros(len(public), dtype=bool)
    metrics: dict[str, object] = {}
    models: dict[str, ExtraTreesRegressor] = {}

    for target_index, target in enumerate(TARGETS):
        labeled = source.dropna(subset=[target])
        positions = np.array(
            [public_position[str(key)] for key in labeled.connectivity_key], dtype=int
        )
        x_source = x_public[positions]
        y_source = labeled[target].to_numpy(dtype=float)
        folds = list(KFold(n_splits=5, shuffle=True, random_state=20260826).split(x_source))
        candidates: list[tuple[float, int, np.ndarray]] = []
        for leaf in (1, 2, 3):
            oof = np.full(len(labeled), np.nan)
            for fold, (train, valid) in enumerate(folds):
                model = estimator(20260826 + 101 * target_index + fold, leaf)
                model.fit(x_source[train], y_source[train])
                oof[valid] = model.predict(x_source[valid])
            candidates.append((float(mean_absolute_error(y_source, oof)), leaf, oof))
        mae, leaf, source_oof = min(candidates, key=lambda item: item[0])
        model = estimator(20260826 + 1009 * target_index, leaf)
        model.fit(x_source, y_source)
        models[target] = model
        public_prediction = model.predict(x_public)
        public_prediction[positions] = source_oof
        public_output[f"{target}_teacher"] = public_prediction
        benchmark_output[f"{target}_teacher"] = model.predict(x_benchmark)
        scope_has_oof[positions] = True
        metrics[target] = {
            "n_source": len(labeled),
            "selected_min_samples_leaf": leaf,
            "oof_mae": mae,
        }
        print(target, metrics[target], flush=True)

    public_output["prediction_scope"] = np.where(
        scope_has_oof, "public_teacher_oof_or_mixed", "public_teacher_full"
    )
    pd.concat([benchmark_output, public_output], ignore_index=True).to_parquet(
        processed / "implicit_solvent_teacher_predictions.parquet", index=False
    )
    model_dir = ROOT / "models/implicit_solvent_teacher"
    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": descriptor_columns,
            "target_columns": list(TARGETS),
        },
        model_dir / "models.joblib",
        compress=3,
    )
    metadata = {
        "source_records_after_benchmark_exclusion": len(source),
        "benchmark_connectivity_overlap": 0,
        "targets": metrics,
        "student_inputs": "RDKit descriptors and Morgan radius-2 fingerprint",
        "inference_simulation": False,
        "source": str(SOURCE_DIR.relative_to(ROOT)),
        "note": "GBn2 and 280k-frame learned-solvent simulation outputs are training-only",
        "model_artifact": str((model_dir / "models.joblib").relative_to(ROOT)),
    }
    (processed / "implicit_solvent_teacher_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
