#!/usr/bin/env python3
"""Run frozen matched endpoint controls for the SolvAI confirmation package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from confirmatory_common import (
    REPEAT_SEEDS,
    SHUFFLE_SEEDS,
    bootstrap_difference,
    endpoint_model,
    fit_predict,
    load_confirmatory_data,
    metric_record,
)


def folds_for(data, split_seed: int | None):
    if split_seed is None:
        folds = data.benchmark.fold_random.to_numpy(dtype=int)
        for fold in sorted(np.unique(folds)):
            yield int(fold), np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
    else:
        splitter = KFold(n_splits=5, shuffle=True, random_state=split_seed)
        for fold, (train, test) in enumerate(splitter.split(data.benchmark)):
            yield fold, train, test


def predict_feature_set(data, benchmark_x, public_x, split_seed: int | None):
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = data.public.delta_g_exp.to_numpy(dtype=float)
    prediction = np.full(len(data.benchmark), np.nan)
    fold_ids = np.full(len(data.benchmark), -1, dtype=int)
    for fold, train, test in folds_for(data, split_seed):
        prediction[test] = fit_predict(
            public_x, public_truth, benchmark_x, truth, train, test
        )
        fold_ids[test] = fold
    if not np.isfinite(prediction).all():
        raise AssertionError("Incomplete endpoint prediction")
    return prediction, fold_ids


def predict_shuffled(data, split_seed: int | None, permutation_seed: int, repeat_index: int):
    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = data.public.delta_g_exp.to_numpy(dtype=float)
    prediction = np.full(len(data.benchmark), np.nan)
    fold_ids = np.full(len(data.benchmark), -1, dtype=int)
    for fold, train, test in folds_for(data, split_seed):
        rng = np.random.default_rng(permutation_seed + 1000 * repeat_index + fold)
        public_response = data.public_responses["full"][rng.permutation(len(data.public))]
        train_response = data.benchmark_responses["full"][train][
            rng.permutation(len(train))
        ]
        test_response = data.benchmark_responses["full"][test][rng.permutation(len(test))]
        public_x = np.column_stack([data.public_structure, public_response])
        train_x = np.column_stack([data.benchmark_structure[train], train_response])
        test_x = np.column_stack([data.benchmark_structure[test], test_response])
        x_fit = np.vstack([public_x, train_x])
        y_fit = np.concatenate([public_truth, truth[train]])
        weights = np.concatenate(
            [np.ones(len(public_truth)), np.full(len(train), 3.0)]
        )
        fold_predictions = []
        for model_seed in (11, 29, 47):
            model = endpoint_model(model_seed)
            model.fit(x_fit, y_fit, extratreesregressor__sample_weight=weights)
            fold_predictions.append(model.predict(test_x))
        prediction[test] = np.mean(fold_predictions, axis=0)
        fold_ids[test] = fold
    return prediction, fold_ids


def rows_for(data, prediction, fold_ids, method, partition, repeat, split_seed, shuffle_seed=None):
    result = data.benchmark[
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            "delta_g_exp",
        ]
    ].copy()
    result = result.rename(columns={"delta_g_exp": "y_true"})
    result["partition"] = partition
    result["repeat"] = repeat
    result["split_seed"] = split_seed
    result["fold"] = fold_ids
    result["method"] = method
    result["shuffle_seed"] = shuffle_seed
    result["y_pred"] = prediction
    result["residual"] = result.y_true - result.y_pred
    result["absolute_error"] = result.residual.abs()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    parser.add_argument(
        "--release-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--mode", choices=("primary", "repeats", "shuffle", "all"), default="all"
    )
    args = parser.parse_args()
    out = args.release_root / "results" / "confirmatory"
    out.mkdir(parents=True, exist_ok=True)
    data = load_confirmatory_data(args.workspace_root)
    prediction_rows = []

    if args.mode in {"primary", "all"}:
        for method, (benchmark_x, public_x) in data.feature_sets.items():
            prediction, fold_ids = predict_feature_set(data, benchmark_x, public_x, None)
            prediction_rows.append(
                rows_for(data, prediction, fold_ids, method, "primary", -1, None)
            )
            print("primary", method, metric_record(data.benchmark.delta_g_exp, prediction), flush=True)

    if args.mode in {"repeats", "all"}:
        for repeat, split_seed in enumerate(REPEAT_SEEDS):
            for method, (benchmark_x, public_x) in data.feature_sets.items():
                prediction, fold_ids = predict_feature_set(
                    data, benchmark_x, public_x, split_seed
                )
                prediction_rows.append(
                    rows_for(
                        data,
                        prediction,
                        fold_ids,
                        method,
                        "repeat",
                        repeat,
                        split_seed,
                    )
                )
                print(
                    "repeat",
                    repeat,
                    method,
                    metric_record(data.benchmark.delta_g_exp, prediction),
                    flush=True,
                )

    if args.mode in {"shuffle", "all"}:
        partitions = [("primary", -1, None, 0)] + [
            ("repeat", repeat, seed, repeat + 1)
            for repeat, seed in enumerate(REPEAT_SEEDS)
        ]
        for partition, repeat, split_seed, repeat_index in partitions:
            for shuffle_seed in SHUFFLE_SEEDS:
                prediction, fold_ids = predict_shuffled(
                    data, split_seed, shuffle_seed, repeat_index
                )
                prediction_rows.append(
                    rows_for(
                        data,
                        prediction,
                        fold_ids,
                        "Z_shuffled_full_priors",
                        partition,
                        repeat,
                        split_seed,
                        shuffle_seed,
                    )
                )
                print(
                    partition,
                    repeat,
                    "shuffle",
                    shuffle_seed,
                    metric_record(data.benchmark.delta_g_exp, prediction),
                    flush=True,
                )

    predictions = pd.concat(prediction_rows, ignore_index=True)
    suffix = args.mode
    prediction_path = out / f"endpoint_{suffix}_predictions.parquet"
    predictions.to_parquet(prediction_path, index=False)

    metric_rows = []
    group_columns = ["partition", "repeat", "split_seed", "method", "shuffle_seed"]
    for keys, group in predictions.groupby(group_columns, dropna=False, sort=False):
        record = dict(zip(group_columns, keys, strict=True))
        record.update(metric_record(group.y_true.to_numpy(), group.y_pred.to_numpy()))
        record["n"] = len(group)
        metric_rows.append(record)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(out / f"endpoint_{suffix}_metrics.csv", index=False)

    comparison_rows = []
    for (partition, repeat, split_seed), group in predictions.groupby(
        ["partition", "repeat", "split_seed"], dropna=False
    ):
        baseline = group[group.method.eq("A_structure_only")]
        if len(baseline) != 85:
            continue
        baseline = baseline.set_index("molecule_id")
        for (method, shuffle_seed), candidate in group.groupby(
            ["method", "shuffle_seed"], dropna=False
        ):
            if method == "A_structure_only":
                continue
            candidate = candidate.set_index("molecule_id").loc[baseline.index]
            comparison_rows.append(
                {
                    "partition": partition,
                    "repeat": repeat,
                    "split_seed": split_seed,
                    "method": method,
                    "shuffle_seed": shuffle_seed,
                    **bootstrap_difference(
                        baseline.y_true.to_numpy(),
                        candidate.y_pred.to_numpy(),
                        baseline.y_pred.to_numpy(),
                    ),
                }
            )
    comparisons = pd.DataFrame(comparison_rows)
    comparisons.to_csv(out / f"endpoint_{suffix}_paired_comparisons.csv", index=False)

    metadata = {
        "mode": args.mode,
        "benchmark_rows": len(data.benchmark),
        "public_rows": len(data.public),
        "structure_features": data.benchmark_structure.shape[1],
        "full_response_features": data.benchmark_responses["full"].shape[1],
        "feature_dimensions": {
            name: int(values[0].shape[1]) for name, values in data.feature_sets.items()
        },
        "model_seeds": [11, 29, 47],
        "repeat_seeds": list(REPEAT_SEEDS),
        "shuffle_seeds": list(SHUFFLE_SEEDS),
        "prediction_file": prediction_path.name,
    }
    (out / f"endpoint_{suffix}_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
