"""Repeated nested confirmation of the frozen structure-only candidate set.

The candidate feature blocks and model settings were frozen after the first
five-fold experiment.  Each repeat creates a new shuffled outer split.  Within
each outer training partition, feature-block selection uses only inner-OOF
errors; the outer-test experimental labels are never consulted.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from arrow_distill.experiments import mae
from run_nested_smd_teacher_confirmation import fit_model, load_problem
from sklearn.model_selection import KFold

DEFAULT_REPEAT_SEEDS = (314159, 271828, 161803, 141421, 173205)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inner-trees", type=int, default=120)
    parser.add_argument("--outer-trees", type=int, default=360)
    parser.add_argument("--model-seeds", type=int, nargs="+", default=[11, 29, 47])
    parser.add_argument("--repeat-seeds", type=int, nargs="+", default=list(DEFAULT_REPEAT_SEEDS))
    parser.add_argument("--output-stem", default="repeated_nested_smd_confsolv_confirmation")
    args = parser.parse_args()

    benchmark, public, feature_sets, source_mask = load_problem(
        no_static_arrow=True,
        confsolv_response=True,
    )
    truth = benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = public.delta_g_exp.to_numpy(dtype=float)[source_mask]
    prediction_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []

    for repeat, split_seed in enumerate(args.repeat_seeds):
        outer_splitter = KFold(n_splits=5, shuffle=True, random_state=split_seed)
        fixed = {name: np.full(len(benchmark), np.nan) for name in feature_sets}
        nested = np.full(len(benchmark), np.nan)
        selected_by_molecule = np.full(len(benchmark), "", dtype=object)
        outer_fold_by_molecule = np.full(len(benchmark), -1, dtype=int)

        for outer_fold, (outer_train, outer_test) in enumerate(outer_splitter.split(benchmark)):
            outer_fold_by_molecule[outer_test] = outer_fold
            inner_splitter = KFold(
                n_splits=4,
                shuffle=True,
                random_state=split_seed + 1000 + outer_fold,
            )
            scores: dict[str, float] = {}
            for name, (benchmark_x, public_x) in feature_sets.items():
                inner_prediction = np.full(len(outer_train), np.nan)
                for inner_fold, (fit_position, valid_position) in enumerate(
                    inner_splitter.split(outer_train)
                ):
                    train = outer_train[fit_position]
                    valid = outer_train[valid_position]
                    seed_predictions = []
                    for model_seed in args.model_seeds:
                        model = fit_model(
                            public_x[source_mask],
                            public_truth,
                            benchmark_x[train],
                            truth[train],
                            model_seed + 100 * inner_fold,
                            args.inner_trees,
                        )
                        seed_predictions.append(model.predict(benchmark_x[valid]))
                    inner_prediction[valid_position] = np.mean(seed_predictions, axis=0)
                scores[name] = mae(truth[outer_train], inner_prediction)

            selected = min(scores, key=scores.get)
            for name, (benchmark_x, public_x) in feature_sets.items():
                fold_predictions = []
                for model_seed in args.model_seeds:
                    model = fit_model(
                        public_x[source_mask],
                        public_truth,
                        benchmark_x[outer_train],
                        truth[outer_train],
                        model_seed,
                        args.outer_trees,
                    )
                    fold_predictions.append(model.predict(benchmark_x[outer_test]))
                fixed[name][outer_test] = np.mean(fold_predictions, axis=0)
            nested[outer_test] = fixed[selected][outer_test]
            selected_by_molecule[outer_test] = selected
            selection_rows.append(
                {
                    "repeat": repeat,
                    "split_seed": split_seed,
                    "outer_fold": outer_fold,
                    "selected": selected,
                    **{f"inner_mae_{name}": score for name, score in scores.items()},
                }
            )
            print(repeat, outer_fold, selected, scores, flush=True)

        candidates = {"Nested selection": nested, **fixed}
        for method, prediction in candidates.items():
            metric_rows.append(
                {
                    "repeat": repeat,
                    "split_seed": split_seed,
                    "method": method,
                    "n": len(benchmark),
                    "mae": mae(truth, prediction),
                }
            )
            for index, value in enumerate(prediction):
                prediction_rows.append(
                    {
                        "repeat": repeat,
                        "split_seed": split_seed,
                        "outer_fold": int(outer_fold_by_molecule[index]),
                        "molecule_id": benchmark.iloc[index].molecule_id,
                        "molecule_name": benchmark.iloc[index].molecule_name,
                        "functional_group_family": benchmark.iloc[index].functional_group_family,
                        "scaffold": benchmark.iloc[index].scaffold,
                        "method": method,
                        "selected_block": (
                            selected_by_molecule[index] if method == "Nested selection" else method
                        ),
                        "y_true": truth[index],
                        "y_pred": float(value),
                        "absolute_error": abs(truth[index] - float(value)),
                    }
                )
        print(
            "repeat",
            repeat,
            {name: mae(truth, values) for name, values in candidates.items()},
            flush=True,
        )

    output_stem = ROOT / "results" / args.output_stem
    pd.DataFrame(prediction_rows).to_parquet(output_stem.with_suffix(".parquet"), index=False)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output_stem.with_suffix(".csv"), index=False)
    pd.DataFrame(selection_rows).to_csv(
        ROOT / "results" / f"{args.output_stem}_selections.csv", index=False
    )
    summary = (
        metrics.groupby("method", as_index=False)
        .agg(
            repeats=("mae", "size"),
            mean_mae=("mae", "mean"),
            std_mae=("mae", "std"),
            min_mae=("mae", "min"),
            max_mae=("mae", "max"),
        )
        .sort_values("mean_mae")
    )
    summary.to_csv(ROOT / "results" / f"{args.output_stem}_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
