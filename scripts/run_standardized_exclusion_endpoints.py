#!/usr/bin/env python3
"""Repeat affected endpoint analyses using standardized-exclusion teachers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from confirmatory_common import (
    MODEL_SEEDS,
    REPEAT_SEEDS,
    bootstrap_difference,
    endpoint_model,
    load_confirmatory_data,
    metric_record,
)
from run_confirmatory_endpoint import predict_feature_set, rows_for


def teacher_overrides(release_root: Path) -> dict[str, Path]:
    root = release_root / "results" / "confirmatory" / "teacher_refits"
    paths = {
        "combisolv_qm": root / "combisolv_qm" / "teacher_predictions.parquet",
        "molsolv_smd": root / "molsolv_smd" / "teacher_predictions.parquet",
        "confsolv": root / "confsolv" / "teacher_predictions.parquet",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing corrected teacher outputs: {missing}")
    return paths


def public_only_predictions(data, method: str) -> np.ndarray:
    benchmark_x, public_x = data.feature_sets[method]
    public_y = data.public.delta_g_exp.to_numpy(dtype=float)
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
    data = load_confirmatory_data(args.workspace_root, teacher_overrides(args.release_root))
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    prediction_rows = []

    for method, (benchmark_x, public_x) in data.feature_sets.items():
        prediction, fold_ids = predict_feature_set(data, benchmark_x, public_x, None)
        prediction_rows.append(
            rows_for(
                data,
                prediction,
                fold_ids,
                method,
                "standardized_exclusion_primary",
                -1,
                None,
            )
        )
        print("primary", method, metric_record(truth, prediction), flush=True)

    for repeat, split_seed in enumerate(REPEAT_SEEDS):
        for method in ("A_structure_only", "F_full_solvai"):
            benchmark_x, public_x = data.feature_sets[method]
            prediction, fold_ids = predict_feature_set(data, benchmark_x, public_x, split_seed)
            prediction_rows.append(
                rows_for(
                    data,
                    prediction,
                    fold_ids,
                    method,
                    "standardized_exclusion_repeat",
                    repeat,
                    split_seed,
                )
            )
            print("repeat", repeat, method, metric_record(truth, prediction), flush=True)

    for method in ("A_structure_only", "F_full_solvai"):
        prediction = public_only_predictions(data, method)
        prediction_rows.append(
            rows_for(
                data,
                prediction,
                np.full(len(data.benchmark), -1, dtype=int),
                method,
                "standardized_exclusion_zero_arrow",
                -1,
                None,
            )
        )
        print("zero_arrow", method, metric_record(truth, prediction), flush=True)

    predictions = pd.concat(prediction_rows, ignore_index=True)
    predictions.to_parquet(out / "standardized_exclusion_endpoint_predictions.parquet", index=False)
    group_columns = ["partition", "repeat", "split_seed", "method"]
    metric_rows = []
    comparison_rows = []
    for keys, group in predictions.groupby(group_columns, dropna=False, sort=False):
        record = dict(zip(group_columns, keys, strict=True))
        record.update(metric_record(group.y_true, group.y_pred))
        record["n"] = len(group)
        metric_rows.append(record)
    for keys, group in predictions.groupby(
        ["partition", "repeat", "split_seed"], dropna=False, sort=False
    ):
        baseline = group.loc[group.method.eq("A_structure_only")].set_index("molecule_id")
        candidate = group.loc[group.method.eq("F_full_solvai")].set_index("molecule_id")
        if len(baseline) != 85 or len(candidate) != 85:
            continue
        candidate = candidate.loc[baseline.index]
        comparison_rows.append(
            {
                "partition": keys[0],
                "repeat": keys[1],
                "split_seed": keys[2],
                **bootstrap_difference(
                    baseline.y_true.to_numpy(),
                    candidate.y_pred.to_numpy(),
                    baseline.y_pred.to_numpy(),
                ),
            }
        )
    pd.DataFrame(metric_rows).to_csv(
        out / "standardized_exclusion_endpoint_metrics.csv", index=False
    )
    pd.DataFrame(comparison_rows).to_csv(
        out / "standardized_exclusion_endpoint_paired_comparisons.csv", index=False
    )

    original_tables = {
        "combisolv_qm": args.workspace_root
        / "data/processed/combisolv_qm_teacher_predictions.parquet",
        "molsolv_smd": args.workspace_root
        / "data/processed/molsolv_smd_teacher_predictions.parquet",
        "confsolv": args.workspace_root
        / "data/processed/confsolv_water_teacher_predictions.parquet",
    }
    shift_rows = []
    for source, corrected_path in teacher_overrides(args.release_root).items():
        old = pd.read_parquet(original_tables[source]).set_index("molecule_id")
        corrected = pd.read_parquet(corrected_path).set_index("molecule_id")
        columns = sorted(
            column for column in old.columns if column.endswith("_teacher") and column in corrected
        )
        for scope, ids in (
            ("benchmark", data.benchmark.molecule_id.astype(str)),
            ("public_endpoint", data.public.molecule_id.astype(str)),
        ):
            for column in columns:
                difference = corrected.reindex(ids)[column].to_numpy(dtype=float) - old.reindex(
                    ids
                )[column].to_numpy(dtype=float)
                shift_rows.append(
                    {
                        "source": source,
                        "scope": scope,
                        "feature": column,
                        "n": len(ids),
                        "mean_change": float(difference.mean()),
                        "mean_absolute_change": float(np.abs(difference).mean()),
                        "maximum_absolute_change": float(np.abs(difference).max()),
                        "correlation": float(
                            np.corrcoef(
                                corrected.reindex(ids)[column].to_numpy(dtype=float),
                                old.reindex(ids)[column].to_numpy(dtype=float),
                            )[0, 1]
                        ),
                    }
                )
    pd.DataFrame(shift_rows).to_csv(
        out / "standardized_exclusion_teacher_prediction_shifts.csv", index=False
    )
    metadata = {
        "teacher_overrides": {
            key: str(value.relative_to(args.release_root))
            for key, value in teacher_overrides(args.release_root).items()
        },
        "model_seeds": list(MODEL_SEEDS),
        "repeat_seeds": list(REPEAT_SEEDS),
        "standardized_exclusions": {
            source: json.loads(
                (
                    args.release_root
                    / "results/confirmatory/teacher_refits"
                    / source
                    / "exclusions.json"
                ).read_text()
            )["count"]
            for source in ("combisolv_qm", "molsolv_smd", "confsolv")
        },
    }
    (out / "standardized_exclusion_endpoint_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
