"""Build leakage-safe structure-only teachers from modern OpenFF ASFE data.

The primary source contains 603 explicit-solvent OpenFF 2.3.0 calculations.
All ARROW benchmark connectivities are excluded before fitting.  For public
experimental rows that also carry a simulation label, teacher predictions are
cross-fitted so the downstream student never receives an in-sample simulation
feature during training.
"""

from __future__ import annotations

import bz2
import glob
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from arrow_distill.data import FREESOLV_JSON, ROOT, canonicalize
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold

CURRENT_RESULTS = (
    ROOT / "data/external/repos/alchemical-benchmark-resources/submissions/"
    "2026_03_17_openff-2.3.0_freesolv/get_results/output/"
    "computational_results.json.bz2"
)
ASH_RESULTS = ROOT / "data/external/repos/ash-sage-rc2/04_benchmark/sfes/results/freesolv"


def magnitude(value: dict[str, object]) -> float:
    if value.get("unit") != "kilocalorie_per_mole":
        raise ValueError(f"Unexpected unit: {value.get('unit')}")
    return float(value["magnitude"])


def load_current(benchmark_keys: set[str]) -> pd.DataFrame:
    freesolv = json.loads(FREESOLV_JSON.read_text())
    with bz2.open(CURRENT_RESULTS, "rt") as stream:
        records = json.load(stream)["dg"]
    rows: list[dict[str, object]] = []
    for record in records:
        public_id = str(record["ligand"])
        source = freesolv[public_id]
        smiles, inchi_key, connectivity_key, _ = canonicalize(source["smiles"])
        if connectivity_key in benchmark_keys:
            continue
        solvent_repeats = np.array([magnitude(item) for item in record["dgs_solvent"]], dtype=float)
        vacuum_repeats = np.array([magnitude(item) for item in record["dgs_vacuum"]], dtype=float)
        dg = magnitude(record["dg"])
        rows.append(
            {
                "public_id": public_id,
                "canonical_smiles": smiles,
                "inchi_key": inchi_key,
                "connectivity_key": connectivity_key,
                "openff23_dg": dg,
                "openff23_uncertainty": magnitude(record["dg_uncertainty"]),
                "openff23_solvent_decoupling": float(solvent_repeats.mean()),
                "openff23_vacuum_decoupling": float(vacuum_repeats.mean()),
                "openff23_repeat_spread": float(np.std(vacuum_repeats - solvent_repeats, ddof=1)),
                "openff23_exp_residual": float(source["expt"]) - dg,
                "openff23_minus_legacy_gaff": dg - float(source["calc"]),
            }
        )
    frame = pd.DataFrame(rows)
    numeric = [
        column
        for column in frame
        if column not in {"public_id", "canonical_smiles", "inchi_key", "connectivity_key"}
    ]
    aggregations: dict[str, object] = {column: "median" for column in numeric}
    aggregations.update(
        {
            "public_id": lambda values: ";".join(sorted(values)),
            "canonical_smiles": "first",
            "inchi_key": "first",
        }
    )
    return frame.groupby("connectivity_key", as_index=False).agg(aggregations)


def load_ash_campaigns(benchmark_keys: set[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for filename in sorted(glob.glob(str(ASH_RESULTS / "*.csv"))):
        path = Path(filename)
        if path.name == "experiment.csv":
            continue
        frame = pd.read_csv(path)
        method = path.stem.replace("-", "_").replace(".", "_")
        for record in frame.itertuples(index=False):
            smiles, inchi_key, connectivity_key, _ = canonicalize(record[3])
            if connectivity_key in benchmark_keys:
                continue
            rows.append(
                {
                    "canonical_smiles": smiles,
                    "inchi_key": inchi_key,
                    "connectivity_key": connectivity_key,
                    "method": method,
                    "value": float(record[5]),
                }
            )
    long = pd.DataFrame(rows)
    pivot = long.pivot_table(
        index="connectivity_key", columns="method", values="value", aggfunc="median"
    ).add_prefix("ash_")
    identity = long.groupby("connectivity_key", as_index=False).agg(
        canonical_smiles=("canonical_smiles", "first"),
        inchi_key=("inchi_key", "first"),
    )
    frame = identity.merge(pivot.reset_index(), on="connectivity_key", how="left")
    method_columns = [column for column in frame if column.startswith("ash_")]
    frame["ash_consensus_dg"] = frame[method_columns].mean(axis=1)
    frame["ash_method_spread"] = frame[method_columns].std(axis=1)
    return frame


def estimator(seed: int, leaf: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=500,
        max_features=0.7,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )


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
    if not benchmark.molecule_id.equals(benchmark_features.molecule_id):
        raise AssertionError("Benchmark feature order changed")
    if not public.molecule_id.equals(public_features.molecule_id):
        raise AssertionError("Public feature order changed")

    benchmark_keys = set(benchmark.inchi_connectivity_key)
    current = load_current(benchmark_keys)
    ash = load_ash_campaigns(benchmark_keys)
    if (set(current.connectivity_key) | set(ash.connectivity_key)) & benchmark_keys:
        raise AssertionError("ARROW benchmark connectivity leaked into OpenFF teachers")
    source = current.merge(
        ash.drop(columns=["canonical_smiles", "inchi_key"]),
        on="connectivity_key",
        how="left",
    )
    source.to_parquet(processed / "openff_alchemical_nonbenchmark.parquet", index=False)

    descriptor_columns = [
        column for column in benchmark_features if column.startswith(("rdkit__", "morgan2__"))
    ]
    x_benchmark = benchmark_features[descriptor_columns].to_numpy(dtype=np.float32)
    x_public = public_features[descriptor_columns].to_numpy(dtype=np.float32)
    public_position = {
        key: index for index, key in enumerate(public.inchi_connectivity_key.astype(str))
    }
    target_columns = [
        "openff23_dg",
        "openff23_solvent_decoupling",
        "openff23_vacuum_decoupling",
        "openff23_exp_residual",
        "openff23_minus_legacy_gaff",
        "ash_consensus_dg",
        "ash_method_spread",
    ]
    benchmark_output = benchmark[
        ["molecule_id", "canonical_smiles", "inchi_connectivity_key"]
    ].rename(columns={"inchi_connectivity_key": "connectivity_key"})
    benchmark_output["prediction_scope"] = "benchmark_external_teacher"
    public_output = public[["molecule_id", "canonical_smiles", "inchi_connectivity_key"]].rename(
        columns={"inchi_connectivity_key": "connectivity_key"}
    )
    scope_has_oof = np.zeros(len(public), dtype=bool)
    metrics: dict[str, object] = {}
    models: dict[str, ExtraTreesRegressor] = {}

    for target_index, target in enumerate(target_columns):
        labeled = source.dropna(subset=[target]).copy()
        positions = np.array(
            [public_position[str(key)] for key in labeled.connectivity_key], dtype=int
        )
        x_source = x_public[positions]
        y_source = labeled[target].to_numpy(dtype=float)
        n_splits = min(5, len(labeled))
        folds = list(KFold(n_splits=n_splits, shuffle=True, random_state=20260826).split(x_source))
        candidates: list[tuple[float, int, np.ndarray]] = []
        for leaf in (1, 2, 3):
            oof = np.full(len(labeled), np.nan)
            for fold, (train, valid) in enumerate(folds):
                model = estimator(20260826 + 101 * target_index + fold, leaf)
                model.fit(x_source[train], y_source[train])
                oof[valid] = model.predict(x_source[valid])
            candidates.append((float(mean_absolute_error(y_source, oof)), leaf, oof))
        mae, selected_leaf, source_oof = min(candidates, key=lambda item: item[0])
        model = estimator(20260826 + 1009 * target_index, selected_leaf)
        model.fit(x_source, y_source)
        models[target] = model
        public_prediction = model.predict(x_public)
        public_prediction[positions] = source_oof
        benchmark_prediction = model.predict(x_benchmark)
        public_output[f"{target}_teacher"] = public_prediction
        benchmark_output[f"{target}_teacher"] = benchmark_prediction
        scope_has_oof[positions] = True
        metrics[target] = {
            "n_source": len(labeled),
            "selected_min_samples_leaf": selected_leaf,
            "oof_mae": mae,
        }
        print(target, metrics[target], flush=True)

    public_output["prediction_scope"] = np.where(
        scope_has_oof, "public_teacher_oof_or_mixed", "public_teacher_full"
    )
    predictions = pd.concat([benchmark_output, public_output], ignore_index=True)
    predictions.to_parquet(processed / "openff_alchemical_teacher_predictions.parquet", index=False)
    model_dir = ROOT / "models/openff_alchemical_teacher"
    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": descriptor_columns,
            "target_columns": target_columns,
        },
        model_dir / "models.joblib",
        compress=3,
    )
    metadata = {
        "openff23_source_records_after_benchmark_exclusion": len(current),
        "ash_source_records_after_benchmark_exclusion": len(ash),
        "benchmark_connectivity_overlap": 0,
        "targets": metrics,
        "student_inputs": "RDKit descriptors and Morgan radius-2 fingerprint",
        "inference_simulation": False,
        "primary_source": str(CURRENT_RESULTS.relative_to(ROOT)),
        "primary_source_license": "CC-BY-4.0",
        "protocol": "OpenFF 2.3.0/NAGL/TIP3P; 14 lambda states; explicit-solvent ASFE",
        "model_artifact": str((model_dir / "models.joblib").relative_to(ROOT)),
    }
    (processed / "openff_alchemical_teacher_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
