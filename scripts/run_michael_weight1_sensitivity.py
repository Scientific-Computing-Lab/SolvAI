#!/usr/bin/env python3
"""Run the prospectively frozen equal-ARROW-weight endpoint sensitivity."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from confirmatory_common import (
    BOOTSTRAP_SEED,
    MODEL_SEEDS,
    bootstrap_difference,
    endpoint_model,
    load_confirmatory_data,
    metric_record,
)
from run_confirmatory_endpoint import folds_for, rows_for
from run_standardized_exclusion_endpoints import teacher_overrides

EXPECTED_HASHES = {
    "data/processed/arrow_solvation_master.parquet": (
        "2b7928f162d094e7ee10d197e66636ba4ae09b0f76d626136b79c0975d3b0310"
    ),
    "data/processed/expanded_public_hydration_nonbenchmark.parquet": (
        "603ed02b6be25d9a3057e321f2c6ea135b012666cfdb8a1b160e37f347951ec4"
    ),
    "data/processed/expanded_public_hydration_rdkit_morgan_features.parquet": (
        "f6d9cd37a90bfc0718261f7251c70100be15947f73d7db12c85659ffc05b28e9"
    ),
    "data/processed/rdkit_morgan_features.parquet": (
        "39877c3938616445d9996e093f8f37744a9c13d56202db9e294466396df298b4"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs(workspace_root: Path) -> None:
    for relative, expected in EXPECTED_HASHES.items():
        observed = sha256(workspace_root / relative)
        if observed != expected:
            raise AssertionError(f"Frozen input hash mismatch: {relative}")


def predict_equal_weight(data, benchmark_x: np.ndarray, public_x: np.ndarray):
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = data.public.delta_g_exp.to_numpy(dtype=float)
    prediction = np.full(len(data.benchmark), np.nan)
    fold_ids = np.full(len(data.benchmark), -1, dtype=int)
    for fold, train, test in folds_for(data, None):
        x_fit = np.vstack([public_x, benchmark_x[train]])
        y_fit = np.concatenate([public_truth, truth[train]])
        weights = np.ones(len(y_fit), dtype=float)
        member_predictions = []
        for seed in MODEL_SEEDS:
            model = endpoint_model(seed)
            model.fit(x_fit, y_fit, extratreesregressor__sample_weight=weights)
            member_predictions.append(model.predict(benchmark_x[test]))
        prediction[test] = np.mean(member_predictions, axis=0)
        fold_ids[test] = fold
    if not np.isfinite(prediction).all():
        raise AssertionError("Incomplete weight-1 prediction")
    return prediction, fold_ids


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    parser.add_argument(
        "--release-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    verify_inputs(args.workspace_root)
    output = args.release_root / "results" / "michael_30aug_sensitivity"
    output.mkdir(parents=True, exist_ok=True)

    data = load_confirmatory_data(
        args.workspace_root, teacher_overrides(args.release_root)
    )
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    prediction_rows = []
    stored = {}
    for method in ("A_structure_only", "F_full_solvai"):
        benchmark_x, public_x = data.feature_sets[method]
        prediction, fold_ids = predict_equal_weight(data, benchmark_x, public_x)
        stored[method] = prediction
        prediction_rows.append(
            rows_for(
                data,
                prediction,
                fold_ids,
                method,
                "michael_30aug_arrow_weight_1",
                -1,
                None,
            )
        )

    predictions = pd.concat(prediction_rows, ignore_index=True)
    predictions.to_parquet(output / "weight1_predictions.parquet", index=False)
    metrics = pd.DataFrame(
        [
            {"method": method, "n": len(truth), **metric_record(truth, prediction)}
            for method, prediction in stored.items()
        ]
    )
    metrics.to_csv(output / "weight1_metrics.csv", index=False)
    comparison = {
        "candidate": "F_full_solvai",
        "baseline": "A_structure_only",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_draws": 100_000,
        **bootstrap_difference(
            truth,
            stored["F_full_solvai"],
            stored["A_structure_only"],
            seed=BOOTSTRAP_SEED,
            draws=100_000,
        ),
    }
    pd.DataFrame([comparison]).to_csv(
        output / "weight1_paired_comparison.csv", index=False
    )
    metadata = {
        "protocol": "release/MICHAEL_30AUG_SENSITIVITY_FREEZE.md",
        "protocol_commit": "9851ed2154be5a8be1f41f3d07975e20fec9a900",
        "benchmark_rows": len(data.benchmark),
        "public_training_rows": len(data.public),
        "external_weight": 1.0,
        "arrow_outer_training_weight": 1.0,
        "fold_column": "fold_random",
        "model_seeds": list(MODEL_SEEDS),
        "input_sha256": EXPECTED_HASHES,
        "teacher_overrides": {
            key: {
                "path": str(path.relative_to(args.release_root)),
                "sha256": sha256(path),
            }
            for key, path in teacher_overrides(args.release_root).items()
        },
        "result": comparison,
    }
    (output / "weight1_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(metrics.to_string(index=False))
    print(pd.DataFrame([comparison]).to_string(index=False))


if __name__ == "__main__":
    main()
