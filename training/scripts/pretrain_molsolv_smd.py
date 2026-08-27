"""Pretrain a structure encoder on leakage-excluded SMD(water) physics.

The source is the deterministic, chemistry-agnostic subset produced by
``prepare_molsolv_smd_teacher.py``.  No ARROW benchmark connectivity is present.
The resulting checkpoint is used only as an initialization for structure-only
downstream prediction.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd
from arrow_distill.data import ROOT


def find_model(directory: Path) -> Path:
    candidates = sorted(directory.rglob("best.pt"))
    if not candidates:
        raise FileNotFoundError(f"No Chemprop model under {directory}")
    return candidates[-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--foundation", default="CHEMELEON")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    processed = ROOT / "data/processed"
    source = pd.read_parquet(processed / "molsolv_smd_water_nonbenchmark.parquet")
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark_keys = set(
        benchmark.loc[benchmark.solvent.eq("water"), "inchi_connectivity_key"].astype(str)
    )
    overlap = set(source.connectivity_key.astype(str)) & benchmark_keys
    if overlap:
        raise AssertionError(f"Benchmark leakage in MolSolv SMD source: {overlap}")
    frame = source[["canonical_smiles", "connectivity_key", "smd_water_dg"]].rename(
        columns={"canonical_smiles": "smiles", "smd_water_dg": "target_smd_water"}
    )
    frame = frame.sample(frac=1, random_state=20260826).reset_index(drop=True)
    n_validation = round(0.05 * len(frame))
    n_test = round(0.05 * len(frame))
    validation = frame.iloc[:n_validation]
    test = frame.iloc[n_validation : n_validation + n_test]
    training = frame.iloc[n_validation + n_test :]

    run_dir = ROOT / "results/molsolv_smd_pretraining"
    input_dir = run_dir / "input"
    model_dir = run_dir / "model"
    input_dir.mkdir(parents=True, exist_ok=True)
    for name, data in (("train", training), ("validation", validation), ("test", test)):
        data.to_csv(input_dir / f"{name}.csv", index=False)

    command = [
        str(ROOT / ".venv/bin/chemprop"),
        "train",
        "-i",
        str(input_dir / "train.csv"),
        str(input_dir / "validation.csv"),
        str(input_dir / "test.csv"),
        "-s",
        "smiles",
        "--target-columns",
        "target_smd_water",
        "--task-type",
        "regression",
        "--epochs",
        str(args.epochs),
        "--patience",
        "6",
        "--batch-size",
        "256",
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
        args.foundation,
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "--num-workers",
        "6",
        "--pytorch-seed",
        "20260826",
        "--output-dir",
        str(model_dir),
        "-q",
    ]
    if args.force or not list(model_dir.rglob("best.pt")):
        with (run_dir / "train.log").open("w") as stream:
            subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)

    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    prediction_frame = pd.concat(
        [
            pd.DataFrame(
                {
                    "molecule_id": benchmark.molecule_id,
                    "canonical_smiles": benchmark.canonical_smiles,
                    "source_split": "benchmark",
                }
            ),
            pd.DataFrame(
                {
                    "molecule_id": public.molecule_id,
                    "canonical_smiles": public.canonical_smiles,
                    "source_split": "expanded_public_hydration",
                }
            ),
        ],
        ignore_index=True,
    )
    prediction_input = input_dir / "benchmark_and_public.csv"
    prediction_frame[["canonical_smiles"]].rename(columns={"canonical_smiles": "smiles"}).to_csv(
        prediction_input, index=False
    )
    raw_prediction_path = run_dir / "benchmark_and_public_raw_predictions.csv"
    subprocess.run(
        [
            str(ROOT / ".venv/bin/chemprop"),
            "predict",
            "-i",
            str(prediction_input),
            "-s",
            "smiles",
            "--model-paths",
            str(find_model(model_dir)),
            "-o",
            str(raw_prediction_path),
            "--accelerator",
            "gpu",
            "--devices",
            "1",
            "-q",
        ],
        check=True,
    )
    raw_predictions = pd.read_csv(raw_prediction_path)
    prediction_columns = [
        column for column in raw_predictions if column.startswith("target_smd_water")
    ]
    if not prediction_columns:
        raise KeyError(f"Missing SMD prediction in {list(raw_predictions)}")
    prediction_frame["molsolv_smd_teacher"] = raw_predictions[prediction_columns[0]].to_numpy(
        dtype=float
    )
    prediction_frame.to_parquet(processed / "molsolv_smd_teacher_predictions.parquet", index=False)

    metadata = {
        "source": "MolSolv M06-2X/6-31G* SMD(water)",
        "rows": len(frame),
        "training_rows": len(training),
        "validation_rows": len(validation),
        "test_rows": len(test),
        "benchmark_connectivity_overlap": 0,
        "foundation": args.foundation,
        "epochs_cap": args.epochs,
        "model_path": str(find_model(model_dir).relative_to(ROOT)),
        "prediction_rows": len(prediction_frame),
        "prediction_output": "data/processed/molsolv_smd_teacher_predictions.parquet",
        "downstream_inference": "structure only",
    }
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
