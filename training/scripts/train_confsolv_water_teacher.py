"""Train structure-only surrogates for the ConfSolv water-response hierarchy.

The source table is globally benchmark-excluded.  The fitted models therefore
turn SMILES-derived descriptors into privileged COSMO-RS/conformer-response
features without exposing any ARROW benchmark label or simulation observable at
inference.
"""

from __future__ import annotations

import argparse
import json

import joblib
import numpy as np
import pandas as pd
from arrow_distill.data import ROOT, rdkit_descriptor_frame
from lightgbm import LGBMRegressor
from sklearn.model_selection import train_test_split

SEED = 20260826
CORE_TARGETS = [
    "confsolv_water_dg_ensemble",
    "confsolv_water_gsolv_gas_lec",
    "confsolv_water_gsolv_solution_lec",
    "confsolv_water_gsolv_gas_weighted",
    "confsolv_water_gsolv_solution_weighted",
    "confsolv_gas_conformer_correction",
    "confsolv_solution_conformer_correction",
    "confsolv_hydration_conformer_correction",
    "confsolv_water_gsolv_std",
    "confsolv_water_response_mean",
    "confsolv_water_response_std",
]


def make_model(seed: int, trees: int) -> LGBMRegressor:
    return LGBMRegressor(
        objective="regression_l1",
        n_estimators=trees,
        learning_rate=0.04,
        num_leaves=31,
        min_child_samples=20,
        subsample=0.8,
        colsample_bytree=0.55,
        reg_lambda=0.1,
        n_jobs=8,
        random_state=seed,
        verbosity=-1,
    )


def aligned_features(frame: pd.DataFrame, descriptor_columns: list[str]) -> np.ndarray:
    values = frame[descriptor_columns].to_numpy(dtype=np.float32)
    values[np.isinf(values)] = np.nan
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", type=int, default=450)
    parser.add_argument("--force-features", action="store_true")
    parser.add_argument(
        "--morgan-bits",
        type=int,
        default=256,
        help="Keep this many highest-variance Morgan bits in addition to RDKit descriptors.",
    )
    args = parser.parse_args()

    processed = ROOT / "data/processed"
    source = pd.read_parquet(processed / "confsolv_water_nonbenchmark.parquet")
    source = source[source.heavy_atom_count.between(1, 18)].copy()
    source = source[source.confsolv_water_dg_ensemble.between(-80.0, 30.0)].reset_index(drop=True)
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    benchmark_keys = set(benchmark.inchi_connectivity_key.astype(str))
    overlap = set(source.connectivity_key.astype(str)) & benchmark_keys
    if overlap:
        raise AssertionError(f"Benchmark connectivity leaked into ConfSolv source: {overlap}")

    feature_path = processed / "confsolv_water_rdkit_morgan_features.parquet"
    if args.force_features or not feature_path.is_file():
        feature_input = source[["canonical_smiles", "inchi_key"]].rename(
            columns={"inchi_key": "molecule_id"}
        )
        source_features = rdkit_descriptor_frame(feature_input)
        source_features.to_parquet(feature_path, index=False)
    else:
        source_features = pd.read_parquet(feature_path)
    source_feature_index = source_features.drop_duplicates("molecule_id").set_index("molecule_id")
    source_features = source_feature_index.reindex(source.inchi_key.astype(str)).reset_index()

    benchmark_features = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    public_features = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    rdkit_columns = [column for column in benchmark_features if column.startswith("rdkit__")]
    morgan_columns = [column for column in benchmark_features if column.startswith("morgan2__")]
    if args.morgan_bits > 0:
        variances = source_features[morgan_columns].var(axis=0).sort_values(ascending=False)
        selected_morgan = variances.head(args.morgan_bits).index.tolist()
    else:
        selected_morgan = []
    descriptor_columns = [*rdkit_columns, *selected_morgan]
    missing_columns = set(descriptor_columns) - set(source_features)
    if missing_columns:
        raise KeyError(f"Missing source descriptors: {sorted(missing_columns)}")
    source_x = aligned_features(source_features, descriptor_columns)
    benchmark_x = aligned_features(benchmark_features, descriptor_columns)
    public_x = aligned_features(public_features, descriptor_columns)

    training_index, validation_index = train_test_split(
        np.arange(len(source)), test_size=0.10, random_state=SEED
    )
    benchmark_output = benchmark[
        ["molecule_id", "canonical_smiles", "inchi_connectivity_key"]
    ].copy()
    benchmark_output["prediction_scope"] = "benchmark_external_teacher"
    public_output = public[["molecule_id", "canonical_smiles", "inchi_connectivity_key"]].copy()
    public_output["prediction_scope"] = "public_external_teacher"
    validation_metrics: dict[str, dict[str, float]] = {}
    models = {}
    for target_index, target in enumerate(CORE_TARGETS):
        y = source[target].to_numpy(dtype=np.float64)
        finite = np.isfinite(y)
        train = training_index[finite[training_index]]
        valid = validation_index[finite[validation_index]]
        validation_model = make_model(SEED + target_index, args.trees)
        validation_model.fit(source_x[train], y[train])
        prediction = validation_model.predict(source_x[valid])
        correlation = float(np.corrcoef(y[valid], prediction)[0, 1])
        validation_metrics[target] = {
            "n_train": len(train),
            "n_validation": len(valid),
            "mae_kcal_mol": float(np.mean(np.abs(y[valid] - prediction))),
            "correlation": correlation,
        }
        model = make_model(SEED + 1000 + target_index, args.trees)
        model.fit(source_x[finite], y[finite])
        models[target] = model
        benchmark_output[f"{target}_teacher"] = model.predict(benchmark_x)
        public_output[f"{target}_teacher"] = model.predict(public_x)
        print(target, validation_metrics[target], flush=True)

    output = pd.concat([benchmark_output, public_output], ignore_index=True)
    output_path = processed / "confsolv_water_teacher_predictions.parquet"
    output.to_parquet(output_path, index=False)
    model_dir = ROOT / "models/confsolv_water_teacher"
    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {"models": models, "descriptor_columns": descriptor_columns},
        model_dir / "lightgbm_models.joblib",
        compress=3,
    )
    metadata = {
        "source": "ConfSolv DOI 10.5281/zenodo.8292520",
        "source_rows_after_chemistry_filter": len(source),
        "benchmark_connectivity_overlap": 0,
        "target_count": len(CORE_TARGETS),
        "target_columns": CORE_TARGETS,
        "validation": validation_metrics,
        "student_inputs": "SMILES-derived RDKit descriptors and Morgan fingerprint only",
        "inference_simulation": False,
        "output": str(output_path.relative_to(ROOT)),
    }
    (model_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
