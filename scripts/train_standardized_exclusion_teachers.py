#!/usr/bin/env python3
"""Refit affected frozen teachers after conservative identity exclusion."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.model_selection import train_test_split

SEED = 20260826
CONFSOLV_TARGETS = [
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


def excluded_smiles(release_root: Path, source: str) -> set[str]:
    exclusions = pd.read_csv(
        release_root / "audits" / "confirmatory" / "standardized_exclusion_records.csv"
    )
    return set(exclusions.loc[exclusions.source.eq(source), "source_smiles"].astype(str))


def prediction_frame(processed: Path) -> pd.DataFrame:
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")].drop_duplicates("molecule_id")
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    return pd.concat(
        [
            benchmark[["molecule_id", "canonical_smiles"]].assign(source_split="benchmark"),
            public[["molecule_id", "canonical_smiles"]].assign(
                source_split="expanded_public_hydration"
            ),
        ],
        ignore_index=True,
    )


def find_model(directory: Path) -> Path:
    candidates = sorted(directory.rglob("best.pt"))
    if not candidates:
        raise FileNotFoundError(f"No Chemprop best.pt under {directory}")
    return candidates[-1]


def train_chemprop(
    workspace: Path,
    run_dir: Path,
    frame: pd.DataFrame,
    target: str,
    exclusions: set[str],
    *,
    epochs: int,
    validation_fraction: float,
    batch_size: int,
    patience: int,
    rdkit_2d: bool,
) -> Path:
    input_dir = run_dir / "input"
    model_dir = run_dir / "model"
    input_dir.mkdir(parents=True, exist_ok=True)
    original_rows = len(frame)
    frame = frame.sample(frac=1, random_state=SEED).reset_index(drop=True)
    n_validation = round(validation_fraction * original_rows)
    n_test = round(validation_fraction * original_rows)
    splits = {
        "validation": frame.iloc[:n_validation],
        "test": frame.iloc[n_validation : n_validation + n_test],
        "train": frame.iloc[n_validation + n_test :],
    }
    splits = {
        name: split.loc[~split.smiles.astype(str).isin(exclusions)]
        for name, split in splits.items()
    }
    for name, split in splits.items():
        split[["smiles", target]].to_csv(input_dir / f"{name}.csv", index=False)
    command = [
        str(workspace / ".venv" / "bin" / "chemprop"),
        "train",
        "-i",
        str(input_dir / "train.csv"),
        str(input_dir / "validation.csv"),
        str(input_dir / "test.csv"),
        "-s",
        "smiles",
        "--target-columns",
        target,
        "--task-type",
        "regression",
        "--epochs",
        str(epochs),
        "--patience",
        str(patience),
        "--batch-size",
        str(batch_size),
        "--ffn-hidden-dim",
        "256",
        "--ffn-num-layers",
        "2",
        "--dropout",
        "0.1",
        "--init-lr",
        "0.00003",
        "--max-lr",
        "0.0001",
        "--final-lr",
        "0.00001",
        "--from-foundation",
        "CHEMELEON",
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "--num-workers",
        "6" if batch_size == 256 else "4",
        "--pytorch-seed",
        str(SEED),
        "--output-dir",
        str(model_dir),
        "-q",
    ]
    if rdkit_2d:
        command.extend(["--molecule-featurizers", "rdkit_2d"])
    with (run_dir / "train.log").open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    metadata = {
        "original_rows": original_rows,
        "rows_after_exclusion": sum(map(len, splits.values())),
        "excluded_rows": original_rows - sum(map(len, splits.values())),
        "training_rows": len(splits["train"]),
        "validation_rows": len(splits["validation"]),
        "test_rows": len(splits["test"]),
        "target": target,
        "epochs_cap": epochs,
        "patience": patience,
        "batch_size": batch_size,
        "foundation": "CHEMELEON",
        "pytorch_seed": SEED,
    }
    (run_dir / "training_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return find_model(model_dir)


def predict_chemprop(
    workspace: Path,
    run_dir: Path,
    model_path: Path,
    output_column: str,
    processed: Path,
    *,
    rdkit_2d: bool,
) -> None:
    frame = prediction_frame(processed)
    prediction_input = run_dir / "input" / "benchmark_and_public.csv"
    frame[["canonical_smiles"]].rename(columns={"canonical_smiles": "smiles"}).to_csv(
        prediction_input, index=False
    )
    raw_path = run_dir / "benchmark_and_public_raw_predictions.csv"
    command = [
        str(workspace / ".venv" / "bin" / "chemprop"),
        "predict",
        "-i",
        str(prediction_input),
        "-s",
        "smiles",
        "--model-paths",
        str(model_path),
        "-o",
        str(raw_path),
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "-q",
    ]
    if rdkit_2d:
        command.extend(["--molecule-featurizers", "rdkit_2d"])
    subprocess.run(command, check=True)
    raw = pd.read_csv(raw_path)
    value_columns = [column for column in raw if column.startswith("target_")]
    if len(value_columns) != 1:
        raise AssertionError(f"Unexpected Chemprop output columns: {list(raw)}")
    frame[output_column] = raw[value_columns[0]].to_numpy(dtype=float)
    frame.to_parquet(run_dir / "teacher_predictions.parquet", index=False)


def run_combisolv(workspace: Path, release: Path, output_root: Path) -> None:
    processed = workspace / "data" / "processed"
    source = pd.read_parquet(processed / "combisolv_qm_water_nonbenchmark.parquet")
    exclusions = excluded_smiles(release, "combisolv_qm")
    source = source.drop_duplicates("solute_connectivity_key")
    frame = pd.DataFrame(
        {
            "smiles": source.solute_canonical_smiles,
            "target_qm_water": source.delta_g_solv_qm,
        }
    )
    run_dir = output_root / "combisolv_qm"
    model = train_chemprop(
        workspace,
        run_dir,
        frame,
        "target_qm_water",
        exclusions,
        epochs=30,
        validation_fraction=0.10,
        batch_size=64,
        patience=7,
        rdkit_2d=True,
    )
    predict_chemprop(
        workspace,
        run_dir,
        model,
        "combisolv_qm_teacher",
        processed,
        rdkit_2d=True,
    )
    (run_dir / "exclusions.json").write_text(
        json.dumps({"count": len(exclusions), "smiles": sorted(exclusions)}, indent=2) + "\n"
    )


def run_molsolv(workspace: Path, release: Path, output_root: Path) -> None:
    processed = workspace / "data" / "processed"
    source = pd.read_parquet(processed / "molsolv_smd_water_nonbenchmark.parquet")
    exclusions = excluded_smiles(release, "molsolv_smd")
    frame = source[["canonical_smiles", "smd_water_dg"]].rename(
        columns={"canonical_smiles": "smiles", "smd_water_dg": "target_smd_water"}
    )
    run_dir = output_root / "molsolv_smd"
    model = train_chemprop(
        workspace,
        run_dir,
        frame,
        "target_smd_water",
        exclusions,
        epochs=25,
        validation_fraction=0.05,
        batch_size=256,
        patience=6,
        rdkit_2d=False,
    )
    predict_chemprop(
        workspace,
        run_dir,
        model,
        "molsolv_smd_teacher",
        processed,
        rdkit_2d=False,
    )
    (run_dir / "exclusions.json").write_text(
        json.dumps({"count": len(exclusions), "smiles": sorted(exclusions)}, indent=2) + "\n"
    )


def confsolv_model(seed: int, trees: int = 450) -> LGBMRegressor:
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


def run_confsolv(workspace: Path, release: Path, output_root: Path) -> None:
    processed = workspace / "data" / "processed"
    source = pd.read_parquet(processed / "confsolv_water_nonbenchmark.parquet")
    source = source.loc[
        source.heavy_atom_count.between(1, 18)
        & source.confsolv_water_dg_ensemble.between(-80.0, 30.0)
    ].reset_index(drop=True)
    exclusions = excluded_smiles(release, "confsolv")
    retained = ~source.canonical_smiles.astype(str).isin(exclusions).to_numpy()
    source_features = (
        pd.read_parquet(processed / "confsolv_water_rdkit_morgan_features.parquet")
        .drop_duplicates("molecule_id")
        .set_index("molecule_id")
    )
    source_features = source_features.reindex(source.inchi_key.astype(str)).reset_index()
    benchmark_features = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    public_features = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    rdkit_columns = [column for column in benchmark_features if column.startswith("rdkit__")]
    morgan_columns = [column for column in benchmark_features if column.startswith("morgan2__")]
    selected_morgan = (
        source_features[morgan_columns]
        .var(axis=0)
        .sort_values(ascending=False)
        .head(256)
        .index.tolist()
    )
    columns = [*rdkit_columns, *selected_morgan]

    def matrix(frame):
        values = frame[columns].to_numpy(dtype=np.float32)
        values[np.isinf(values)] = np.nan
        return values

    source_x = matrix(source_features)
    benchmark_x = matrix(benchmark_features)
    public_x = matrix(public_features)
    train, valid = train_test_split(np.arange(len(source)), test_size=0.10, random_state=SEED)
    output = prediction_frame(processed)
    output["inchi_connectivity_key"] = pd.concat(
        [
            pd.read_parquet(processed / "arrow_solvation_master.parquet")
            .loc[lambda frame: frame.solvent.eq("water")]
            .drop_duplicates("molecule_id")
            .inchi_connectivity_key,
            pd.read_parquet(
                processed / "expanded_public_hydration_nonbenchmark.parquet"
            ).inchi_connectivity_key,
        ],
        ignore_index=True,
    )
    validation = {}
    models = {}
    for target_index, target in enumerate(CONFSOLV_TARGETS):
        y = source[target].to_numpy(dtype=float)
        finite = np.isfinite(y) & retained
        target_train = train[finite[train]]
        target_valid = valid[finite[valid]]
        check = confsolv_model(SEED + target_index)
        check.fit(source_x[target_train], y[target_train])
        prediction = check.predict(source_x[target_valid])
        validation[target] = {
            "n_train": len(target_train),
            "n_validation": len(target_valid),
            "mae_kcal_mol": float(np.mean(np.abs(y[target_valid] - prediction))),
            "correlation": float(np.corrcoef(y[target_valid], prediction)[0, 1]),
        }
        model = confsolv_model(SEED + 1000 + target_index)
        model.fit(source_x[finite], y[finite])
        models[target] = model
        values = np.concatenate([model.predict(benchmark_x), model.predict(public_x)])
        output[f"{target}_teacher"] = values
        print(target, validation[target], flush=True)
    run_dir = output_root / "confsolv"
    run_dir.mkdir(parents=True, exist_ok=True)
    output.to_parquet(run_dir / "teacher_predictions.parquet", index=False)
    joblib.dump(
        {"models": models, "descriptor_columns": columns},
        run_dir / "lightgbm_models.joblib",
        compress=3,
    )
    (run_dir / "training_metadata.json").write_text(
        json.dumps(
            {
                "original_rows": len(source),
                "rows_after_exclusion": int(retained.sum()),
                "exclusions": len(exclusions),
                "targets": CONFSOLV_TARGETS,
                "validation": validation,
                "seed": SEED,
                "trees": 450,
                "morgan_bits": 256,
            },
            indent=2,
        )
        + "\n"
    )
    (run_dir / "exclusions.json").write_text(
        json.dumps({"count": len(exclusions), "smiles": sorted(exclusions)}, indent=2) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--release-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--stage", choices=("combisolv", "molsolv", "confsolv", "all"), default="all"
    )
    args = parser.parse_args()
    output_root = args.release_root / "results" / "confirmatory" / "teacher_refits"
    output_root.mkdir(parents=True, exist_ok=True)
    if args.stage in {"combisolv", "all"}:
        run_combisolv(args.workspace_root, args.release_root, output_root)
    if args.stage in {"molsolv", "all"}:
        run_molsolv(args.workspace_root, args.release_root, output_root)
    if args.stage in {"confsolv", "all"}:
        run_confsolv(args.workspace_root, args.release_root, output_root)


if __name__ == "__main__":
    main()
