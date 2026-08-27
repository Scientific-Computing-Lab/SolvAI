"""Ingest the SoluteML experimental hydration and Abraham-parameter archive.

The source uses a 1 M gas / 1 M liquid standard state at 298 K.  Every row
sharing an InChI connectivity block with the ARROW benchmark is globally removed.
Abraham E/S/A/B/L values are used only as training-time targets: downstream
benchmark and public rows receive cross-fitted or out-of-source predictions from
molecular structure.
"""

from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from arrow_distill.data import ROOT, canonicalize, rdkit_descriptor_frame
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold

SOURCE = ROOT / "data/external/zenodo/5792296/Solvation_data-1.0.0/selected_data_for_this_work"
ABRAHAM_TARGETS = ("abraham_e", "abraham_s", "abraham_a", "abraham_b", "abraham_l")


def canonical_rows(frame: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
    identities = []
    for smiles in frame[smiles_column].astype(str):
        canonical, inchi_key, connectivity_key, stereo = canonicalize(smiles)
        identities.append((canonical, inchi_key, connectivity_key, stereo))
    result = frame.copy()
    result[["canonical_smiles", "inchi_key", "connectivity_key", "has_stereo"]] = identities
    return result


def load_hydration(benchmark_keys: set[str]) -> pd.DataFrame:
    path = SOURCE / "dGsolvDB3_selected_data.xlsx"
    frame = pd.read_excel(path, sheet_name="data")
    frame = frame[frame.smiles_solvent.astype(str).eq("O")].copy()
    frame = canonical_rows(frame, "smiles_solute")
    frame = frame[~frame.connectivity_key.isin(benchmark_keys)].copy()
    frame = frame.rename(
        columns={
            "dGsolv_avg [kcal/mol]": "delta_g_exp",
            "dGsolv_std [kcal/mol]": "delta_g_exp_uncertainty",
            "no. of data": "source_measurement_count",
            "Source_all": "source_references",
        }
    )
    frame["molecule_id"] = frame.inchi_key
    frame["source_dataset"] = "SoluteML dGsolvDB3"
    aggregations: dict[str, object] = {
        "molecule_id": "first",
        "canonical_smiles": "first",
        "inchi_key": "first",
        "has_stereo": "max",
        "delta_g_exp": "median",
        "delta_g_exp_uncertainty": "median",
        "source_measurement_count": "sum",
        "source_references": lambda values: "; ".join(sorted(set(map(str, values)))),
        "source_dataset": "first",
    }
    return frame.groupby("connectivity_key", as_index=False).agg(aggregations)


def load_abraham(benchmark_keys: set[str]) -> pd.DataFrame:
    frame = pd.read_excel(SOURCE / "SoluteDB_selected_data.xlsx")
    frame = canonical_rows(frame, "SMILES")
    frame = frame[~frame.connectivity_key.isin(benchmark_keys)].copy()
    frame = frame.rename(
        columns={
            "E": "abraham_e",
            "S": "abraham_s",
            "A": "abraham_a",
            "B": "abraham_b",
            "L": "abraham_l",
        }
    )
    frame["molecule_id"] = frame.inchi_key
    aggregations: dict[str, object] = {
        "molecule_id": "first",
        "canonical_smiles": "first",
        "inchi_key": "first",
        "has_stereo": "max",
    }
    aggregations.update({target: "median" for target in ABRAHAM_TARGETS})
    return frame.groupby("connectivity_key", as_index=False).agg(aggregations)


def estimator(seed: int, leaf: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=160,
        max_features=0.7,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )


def main() -> None:
    processed = ROOT / "data/processed"
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    current_public = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
    benchmark_keys = set(benchmark.inchi_connectivity_key.astype(str))

    hydration = load_hydration(benchmark_keys)
    abraham = load_abraham(benchmark_keys)
    if (set(hydration.connectivity_key) | set(abraham.connectivity_key)) & benchmark_keys:
        raise AssertionError("Benchmark connectivity leaked into SoluteML artifacts")
    hydration.to_parquet(processed / "soluteml_hydration_nonbenchmark.parquet", index=False)
    abraham.to_parquet(processed / "soluteml_abraham_nonbenchmark.parquet", index=False)

    current = current_public.copy()
    current["source_dataset"] = current.get("source_dataset", "existing_public")
    additions = hydration[~hydration.connectivity_key.isin(current.inchi_connectivity_key)].copy()
    additions["inchi_connectivity_key"] = additions.connectivity_key
    common_columns = sorted(set(current) | set(additions))
    expanded = pd.concat(
        [current.reindex(columns=common_columns), additions.reindex(columns=common_columns)],
        ignore_index=True,
    )
    if expanded.inchi_connectivity_key.duplicated().any():
        raise AssertionError("Expanded hydration set contains duplicate connectivities")
    expanded.to_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet", index=False)
    expanded_feature_path = processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    if expanded_feature_path.exists():
        expanded_features = pd.read_parquet(expanded_feature_path)
        if not expanded.molecule_id.equals(expanded_features.molecule_id):
            raise AssertionError("Cached expanded feature order changed")
    else:
        expanded_features = rdkit_descriptor_frame(expanded)
        expanded_features.to_parquet(expanded_feature_path, index=False)

    destination = pd.concat(
        [
            benchmark[["molecule_id", "canonical_smiles", "inchi_connectivity_key"]],
            expanded[["molecule_id", "canonical_smiles", "inchi_connectivity_key"]],
        ],
        ignore_index=True,
    )
    destination_features = pd.concat(
        [
            pd.read_parquet(processed / "rdkit_morgan_features.parquet"),
            expanded_features,
        ],
        ignore_index=True,
    )
    source_feature_path = processed / "soluteml_abraham_rdkit_morgan_features.parquet"
    if source_feature_path.exists():
        source_features = pd.read_parquet(source_feature_path)
        if not abraham.molecule_id.equals(source_features.molecule_id):
            raise AssertionError("Cached Abraham feature order changed")
    else:
        source_features = rdkit_descriptor_frame(abraham)
        source_features.to_parquet(source_feature_path, index=False)
    feature_columns = [
        column for column in source_features if column.startswith(("rdkit__", "morgan2__"))
    ]
    x_source = source_features[feature_columns].to_numpy(dtype=np.float32)
    x_destination = destination_features[feature_columns].to_numpy(dtype=np.float32)
    destination_positions = {
        key: index for index, key in enumerate(destination.inchi_connectivity_key.astype(str))
    }
    output = destination.rename(columns={"inchi_connectivity_key": "connectivity_key"}).copy()
    output["prediction_scope"] = "external_teacher_full"
    metrics: dict[str, object] = {}
    models: dict[str, ExtraTreesRegressor] = {}
    for target_index, target in enumerate(ABRAHAM_TARGETS):
        labeled = np.flatnonzero(abraham[target].notna().to_numpy())
        x_labeled = x_source[labeled]
        y_labeled = abraham[target].iloc[labeled].to_numpy(dtype=float)
        folds = list(KFold(5, shuffle=True, random_state=20260826).split(x_labeled))
        selected_leaf = 2
        source_oof = np.full(len(labeled), np.nan)
        for fold, (train, valid) in enumerate(folds):
            model = estimator(20260826 + 101 * target_index + fold, selected_leaf)
            model.fit(x_labeled[train], y_labeled[train])
            source_oof[valid] = model.predict(x_labeled[valid])
        mae = mean_absolute_error(y_labeled, source_oof)
        model = estimator(20260826 + 1009 * target_index, selected_leaf)
        model.fit(x_labeled, y_labeled)
        models[target] = model
        prediction = model.predict(x_destination)
        for source_index, value in zip(labeled, source_oof, strict=True):
            key = str(abraham.connectivity_key.iloc[source_index])
            if key in destination_positions:
                position = destination_positions[key]
                prediction[position] = value
                output.loc[position, "prediction_scope"] = "external_teacher_oof_or_mixed"
        output[f"{target}_teacher"] = prediction
        metrics[target] = {
            "n_source": len(labeled),
            "selected_min_samples_leaf": selected_leaf,
            "oof_mae": float(mae),
        }
        print(target, metrics[target], flush=True)
    output.to_parquet(processed / "soluteml_abraham_teacher_predictions.parquet", index=False)
    model_dir = ROOT / "models/soluteml_abraham_teacher"
    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": feature_columns,
            "target_columns": list(ABRAHAM_TARGETS),
        },
        model_dir / "models.joblib",
        compress=3,
    )

    metadata = {
        "source": "Zenodo 5792296 Solvation_data-1.0.0",
        "hydration_rows_after_benchmark_exclusion": len(hydration),
        "expanded_hydration_rows": len(expanded),
        "new_hydration_connectivities": len(additions),
        "abraham_rows_after_benchmark_exclusion": len(abraham),
        "benchmark_connectivity_overlap": 0,
        "abraham_targets": metrics,
        "standard_state": "1 mol/L gas to 1 mol/L liquid at 298 K",
        "inference_simulation": False,
        "model_artifact": str((model_dir / "models.joblib").relative_to(ROOT)),
    }
    (processed / "soluteml_hydration_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
