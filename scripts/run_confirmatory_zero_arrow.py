#!/usr/bin/env python3
"""Run the frozen zero-ARROW-label endpoint transfer test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from confirmatory_common import (
    MODEL_SEEDS,
    bootstrap_difference,
    endpoint_model,
    load_confirmatory_data,
    metric_record,
)


def fit_public_only(public_x, public_y, benchmark_x):
    predictions = []
    for seed in MODEL_SEEDS:
        model = endpoint_model(seed)
        model.fit(
            public_x,
            public_y,
            extratreesregressor__sample_weight=np.ones(len(public_y)),
        )
        predictions.append(model.predict(benchmark_x))
    return np.mean(predictions, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--release-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    out = args.release_root / "results" / "confirmatory"
    out.mkdir(parents=True, exist_ok=True)
    data = load_confirmatory_data(args.workspace_root)
    public_truth = data.public.delta_g_exp.to_numpy(dtype=float)
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    rows = []
    stored = {}
    for method in ("A_structure_only", "F_full_solvai"):
        benchmark_x, public_x = data.feature_sets[method]
        prediction = fit_public_only(public_x, public_truth, benchmark_x)
        stored[method] = prediction
        frame = data.benchmark[["molecule_id", "molecule_name", "canonical_smiles"]].copy()
        frame["method"] = method
        frame["y_true"] = truth
        frame["y_pred"] = prediction
        frame["residual"] = truth - prediction
        frame["absolute_error"] = np.abs(frame.residual)
        rows.append(frame)
        print(method, metric_record(truth, prediction), flush=True)
    predictions = pd.concat(rows, ignore_index=True)
    predictions.to_parquet(out / "zero_arrow_transfer_predictions.parquet", index=False)
    metrics = pd.DataFrame(
        [
            {"method": method, "n": len(truth), **metric_record(truth, prediction)}
            for method, prediction in stored.items()
        ]
    )
    metrics.to_csv(out / "zero_arrow_transfer_metrics.csv", index=False)
    comparison = {
        "candidate": "F_full_solvai",
        "baseline": "A_structure_only",
        **bootstrap_difference(truth, stored["F_full_solvai"], stored["A_structure_only"]),
    }
    pd.DataFrame([comparison]).to_csv(
        out / "zero_arrow_transfer_paired_comparison.csv", index=False
    )
    metadata = {
        "benchmark_rows": len(data.benchmark),
        "public_training_rows": len(data.public),
        "arrow_training_labels": 0,
        "model_seeds": list(MODEL_SEEDS),
    }
    (out / "zero_arrow_transfer_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
