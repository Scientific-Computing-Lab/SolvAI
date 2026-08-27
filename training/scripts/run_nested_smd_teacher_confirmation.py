"""Nested confirmation of the narrow SMD-water physics teacher.

This test compares one fixed response model with and without one new source of
information: a structure-predicted SMD(water) free energy trained after global
ARROW-benchmark exclusion.  Feature-block selection occurs only on inner OOF
predictions from the outer training molecules.  The evaluated molecule supplies
SMILES/static descriptors only.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from arrow_distill.experiments import (
    REGIMES,
    extra_trees,
    fit_with_optional_weights,
    inner_splits,
    mae,
    make_prediction_rows,
)


def align(table: pd.DataFrame, ids: pd.Series, columns: list[str]) -> np.ndarray:
    indexed = table.drop_duplicates("molecule_id").set_index("molecule_id")
    values = indexed.reindex(ids.astype(str))[columns].to_numpy(dtype=np.float32)
    if not np.isfinite(values).all():
        raise AssertionError(f"Non-finite aligned values in {columns}")
    return values


def fit_model(
    public_x: np.ndarray,
    public_y: np.ndarray,
    benchmark_x: np.ndarray,
    benchmark_y: np.ndarray,
    seed: int,
    trees: int,
):
    model = extra_trees(
        seed=seed,
        n_estimators=trees,
        max_features=0.7,
        min_samples_leaf=2,
    )
    weights = np.concatenate([np.ones(len(public_y)), np.full(len(benchmark_y), 3.0)])
    return fit_with_optional_weights(
        model,
        np.vstack([public_x, benchmark_x]),
        np.concatenate([public_y, benchmark_y]),
        weights,
    )


def load_problem(*, no_static_arrow: bool, confsolv_response: bool):
    processed = ROOT / "data/processed"
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    features = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    public_features = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    static = pd.read_parquet(processed / "arrow_static_features.parquet")
    qm = pd.read_parquet(processed / "combisolv_qm_teacher_predictions.parquet")
    abraham = pd.read_parquet(processed / "soluteml_abraham_teacher_predictions.parquet")
    openff = pd.read_parquet(processed / "openff_alchemical_teacher_predictions.parquet")
    implicit = pd.read_parquet(processed / "implicit_solvent_teacher_predictions.parquet")
    smd = pd.read_parquet(processed / "molsolv_smd_teacher_predictions.parquet")
    confsolv_path = processed / "confsolv_water_teacher_predictions.parquet"
    confsolv = pd.read_parquet(confsolv_path) if confsolv_response else None

    benchmark_keys = set(benchmark.inchi_connectivity_key.astype(str))
    if benchmark_keys & set(public.inchi_connectivity_key.astype(str)):
        raise AssertionError("Benchmark connectivity leaked into public training")
    descriptor_columns = [
        column for column in features if column.startswith(("rdkit__", "morgan2__"))
    ]
    arrow_columns = [
        column
        for column in static
        if column.startswith("arrow_")
        and column not in {"arrow_hin_path", "arrow_static_available"}
        and pd.api.types.is_numeric_dtype(static[column])
    ]
    if no_static_arrow:
        benchmark_structure = features[descriptor_columns].to_numpy(dtype=np.float32)
        public_structure = public_features[descriptor_columns].to_numpy(dtype=np.float32)
    else:
        benchmark_structure = np.column_stack(
            [
                features[descriptor_columns].to_numpy(dtype=np.float32),
                static.arrow_static_available.astype(float).to_numpy(),
                static[arrow_columns].to_numpy(dtype=np.float32),
            ]
        )
        public_structure = np.column_stack(
            [
                public_features[descriptor_columns].to_numpy(dtype=np.float32),
                np.zeros((len(public), 1 + len(arrow_columns)), dtype=np.float32),
            ]
        )
    qcols = ["combisolv_qm_teacher"]
    acols = [f"abraham_{name}_teacher" for name in "esabl"]
    ocols = ["openff23_dg_teacher", "openff23_exp_residual_teacher"]
    icols = ["gbn2_alchemical_dg_teacher", "gbn2_exp_residual_teacher"]
    benchmark_base = np.column_stack(
        [
            benchmark_structure,
            align(qm, benchmark.molecule_id, qcols),
            align(abraham, benchmark.molecule_id, acols),
            align(openff, benchmark.molecule_id, ocols).sum(axis=1),
            align(implicit, benchmark.molecule_id, icols).sum(axis=1),
        ]
    )
    public_base = np.column_stack(
        [
            public_structure,
            align(qm, public.molecule_id, qcols),
            align(abraham, public.molecule_id, acols),
            align(openff, public.molecule_id, ocols).sum(axis=1),
            align(implicit, public.molecule_id, icols).sum(axis=1),
        ]
    )
    feature_sets = {
        "narrow response without SMD": (benchmark_base, public_base),
        "narrow response + SMD water": (
            np.column_stack(
                [benchmark_base, align(smd, benchmark.molecule_id, ["molsolv_smd_teacher"])]
            ),
            np.column_stack([public_base, align(smd, public.molecule_id, ["molsolv_smd_teacher"])]),
        ),
    }
    if confsolv_response:
        confsolv_columns = [
            "confsolv_gas_conformer_correction_teacher",
            "confsolv_solution_conformer_correction_teacher",
            "confsolv_hydration_conformer_correction_teacher",
            "confsolv_water_gsolv_std_teacher",
            "confsolv_water_response_mean_teacher",
            "confsolv_water_response_std_teacher",
        ]
        feature_sets["narrow response + SMD + ConfSolv response"] = (
            np.column_stack(
                [
                    feature_sets["narrow response + SMD water"][0],
                    align(confsolv, benchmark.molecule_id, confsolv_columns),
                ]
            ),
            np.column_stack(
                [
                    feature_sets["narrow response + SMD water"][1],
                    align(confsolv, public.molecule_id, confsolv_columns),
                ]
            ),
        )
    old_public = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
    old_keys = set(old_public.inchi_connectivity_key.astype(str))
    source_mask = (
        public.inchi_connectivity_key.astype(str).isin(old_keys)
        | public.source_measurement_count.fillna(0).ge(2)
    ).to_numpy()
    return benchmark, public, feature_sets, source_mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regimes", nargs="+", choices=tuple(REGIMES), default=["random_oof"])
    parser.add_argument("--inner-trees", type=int, default=120)
    parser.add_argument("--outer-trees", type=int, default=360)
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 29, 47])
    parser.add_argument("--output-stem", default="nested_smd_teacher_confirmation")
    parser.add_argument("--no-static-arrow", action="store_true")
    parser.add_argument("--confsolv-response", action="store_true")
    args = parser.parse_args()

    benchmark, public, feature_sets, source_mask = load_problem(
        no_static_arrow=args.no_static_arrow,
        confsolv_response=args.confsolv_response,
    )
    truth = benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = public.delta_g_exp.to_numpy(dtype=float)[source_mask]
    rows = []
    selection_rows = []
    for regime in args.regimes:
        folds = benchmark[REGIMES[regime]].to_numpy()
        fixed = {name: np.full(len(benchmark), np.nan) for name in feature_sets}
        nested = np.full(len(benchmark), np.nan)
        details = [""] * len(benchmark)
        for outer_fold in sorted(np.unique(folds)):
            outer_test = np.flatnonzero(folds == outer_fold)
            outer_train = np.flatnonzero(folds != outer_fold)
            inner = inner_splits(outer_train, benchmark, regime)
            scores = {}
            for name, (benchmark_x, public_x) in feature_sets.items():
                prediction = np.full(len(benchmark), np.nan)
                for inner_fold, (inner_train, inner_valid) in enumerate(inner):
                    seed_predictions = []
                    for seed in args.seeds:
                        model = fit_model(
                            public_x[source_mask],
                            public_truth,
                            benchmark_x[inner_train],
                            truth[inner_train],
                            seed + 100 * inner_fold,
                            args.inner_trees,
                        )
                        seed_predictions.append(model.predict(benchmark_x[inner_valid]))
                    prediction[inner_valid] = np.mean(seed_predictions, axis=0)
                valid = np.concatenate([part for _, part in inner])
                scores[name] = mae(truth[valid], prediction[valid])
            selected = min(scores, key=scores.get)
            for name, (benchmark_x, public_x) in feature_sets.items():
                seed_predictions = []
                for seed in args.seeds:
                    model = fit_model(
                        public_x[source_mask],
                        public_truth,
                        benchmark_x[outer_train],
                        truth[outer_train],
                        seed,
                        args.outer_trees,
                    )
                    seed_predictions.append(model.predict(benchmark_x[outer_test]))
                fixed[name][outer_test] = np.mean(seed_predictions, axis=0)
            nested[outer_test] = fixed[selected][outer_test]
            detail = f"selected={selected};" + ";".join(
                f"inner_{name}={score:.6f}" for name, score in sorted(scores.items())
            )
            for index in outer_test:
                details[index] = detail
            selection_rows.append(
                {
                    "regime": regime,
                    "outer_fold": int(outer_fold),
                    "selected": selected,
                    **{f"inner_mae_{name}": score for name, score in scores.items()},
                }
            )
            print(regime, outer_fold, selected, scores, flush=True)

        rows.append(
            make_prediction_rows(
                benchmark,
                nested,
                "Nested narrow SMD teacher selection",
                regime,
                False,
                details=details,
            )
        )
        for name, prediction in fixed.items():
            rows.append(
                make_prediction_rows(
                    benchmark,
                    prediction,
                    f"Fixed {name}",
                    regime,
                    False,
                    details=["structure-only fixed matched ablation"] * len(benchmark),
                )
            )

    output = pd.concat(rows, ignore_index=True)
    path = ROOT / f"results/{args.output_stem}.parquet"
    output.to_parquet(path, index=False)
    metrics = (
        output.groupby(["method", "regime"], as_index=False)
        .agg(n=("absolute_error", "size"), mae=("absolute_error", "mean"))
        .sort_values(["regime", "mae"])
    )
    metrics.to_csv(path.with_suffix(".csv"), index=False)
    pd.DataFrame(selection_rows).to_csv(
        ROOT / f"results/{args.output_stem}_selections.csv", index=False
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
