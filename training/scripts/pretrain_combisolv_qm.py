"""Pretrain a structure encoder on benchmark-excluded CombiSolv-QM water values."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd
from arrow_distill.data import ROOT


def find_model(directory: Path) -> Path:
    candidates = sorted(directory.rglob("*.pt"))
    if not candidates:
        raise FileNotFoundError(f"No Chemprop model under {directory}")
    return candidates[-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--foundation", default="CHEMELEON")
    args = parser.parse_args()
    source = pd.read_parquet(ROOT / "data/processed/combisolv_qm_water_nonbenchmark.parquet")
    source = source.drop_duplicates("solute_connectivity_key").reset_index(drop=True)
    frame = pd.DataFrame(
        {
            "smiles": source.solute_canonical_smiles,
            "target_qm_water": source.delta_g_solv_qm,
        }
    ).sample(frac=1, random_state=20260826)
    n_validation = round(0.1 * len(frame))
    n_test = round(0.1 * len(frame))
    validation = frame.iloc[:n_validation]
    test = frame.iloc[n_validation : n_validation + n_test]
    training = frame.iloc[n_validation + n_test :]
    run_dir = ROOT / "results/combisolv_qm_pretraining"
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
        "target_qm_water",
        "--task-type",
        "regression",
        "--epochs",
        str(args.epochs),
        "--patience",
        "7",
        "--batch-size",
        "64",
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
        "--molecule-featurizers",
        "rdkit_2d",
        "--from-foundation",
        args.foundation,
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "--num-workers",
        "4",
        "--pytorch-seed",
        "20260826",
        "--output-dir",
        str(model_dir),
        "-q",
    ]
    if not list(model_dir.rglob("*.pt")):
        with (run_dir / "train.log").open("w") as stream:
            subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    benchmark = pd.read_parquet(ROOT / "data/processed/arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    expanded_public = ROOT / "data/processed/expanded_public_hydration_nonbenchmark.parquet"
    public = pd.read_parquet(
        expanded_public
        if expanded_public.exists()
        else ROOT / "data/processed/public_hydration_nonbenchmark.parquet"
    )
    prediction_input = input_dir / "benchmark.csv"
    prediction_frame = pd.concat(
        [
            pd.DataFrame(
                {
                    "molecule_id": benchmark.molecule_id,
                    "smiles": benchmark.canonical_smiles,
                    "source_split": "benchmark",
                }
            ),
            pd.DataFrame(
                {
                    "molecule_id": public.molecule_id,
                    "smiles": public.canonical_smiles,
                    "source_split": "expanded_public_hydration",
                }
            ),
        ],
        ignore_index=True,
    )
    prediction_frame[["smiles"]].to_csv(prediction_input, index=False)
    prediction_path = run_dir / "benchmark_predictions.csv"
    predict_command = [
        str(ROOT / ".venv/bin/chemprop"),
        "predict",
        "-i",
        str(prediction_input),
        "-s",
        "smiles",
        "--model-paths",
        str(find_model(model_dir)),
        "--molecule-featurizers",
        "rdkit_2d",
        "-o",
        str(prediction_path),
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "-q",
    ]
    subprocess.run(predict_command, check=True)
    predicted = pd.read_csv(prediction_path)
    value_columns = [column for column in predicted if column.startswith("target_qm_water")]
    output = prediction_frame.rename(columns={"smiles": "canonical_smiles"})
    output["combisolv_qm_teacher"] = predicted[value_columns[0]].to_numpy()
    output.to_parquet(ROOT / "data/processed/combisolv_qm_teacher_predictions.parquet", index=False)
    metadata = {
        "training_rows": len(training),
        "validation_rows": len(validation),
        "test_rows": len(test),
        "benchmark_overlap": int(source.benchmark_solute_overlap.sum()),
        "foundation": args.foundation,
        "epochs_cap": args.epochs,
        "model_path": str(find_model(model_dir).relative_to(ROOT)),
    }
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
