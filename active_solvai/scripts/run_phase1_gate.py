#!/usr/bin/env python3
"""Execute the preregistered Active SolvAI actual-observation gate."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error

RELEASE_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = RELEASE_ROOT.parents[1]
sys.path.insert(0, str(RELEASE_ROOT / "scripts"))

from confirmatory_common import endpoint_model, load_confirmatory_data
from run_standardized_exclusion_endpoints import teacher_overrides

from active_solvai.gate import (
    MODEL_SEEDS,
    apply_affine,
    choose_ridge_alpha,
    condition_gaussian,
    ensemble_fit_predict,
    fit_affine,
    fit_gaussian_curve,
    integration_proxy,
    ridge_correct,
)
from active_solvai.ledger import append_record, sha256

BOOTSTRAP_SEED = 20260828
SHUFFLE_SEEDS = (88001, 88002, 88003, 88004, 88005)
TRAJECTORY_FRACTIONS = (0.1, 0.2, 0.4, 0.7, 1.0)
TOTAL_COMPONENT = "dHdL"
COMPONENTS = (
    "dHdL_Coul",
    "dHdL_Coul_SC",
    "dHdL_VdW",
    "dHdL_VdW_SC",
    "dHdL_LRCor",
    "dHdL_PME",
)
LAMBDA_LABELS = ("0.1", "0.5", "0.9")
SUBSETS = tuple(subset for size in (1, 2, 3) for subset in combinations(range(3), size))


@dataclass
class ResponseArrays:
    total: np.ndarray
    components: np.ndarray
    total_sem: np.ndarray


def key_for_subset(subset: tuple[int, ...]) -> str:
    return "_".join(LAMBDA_LABELS[index].replace(".", "p") for index in subset)


def metric_record(truth: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = np.asarray(truth) - np.asarray(prediction)
    return {
        "mae": float(mean_absolute_error(truth, prediction)),
        "rmse": float(mean_squared_error(truth, prediction) ** 0.5),
        "median_absolute_error": float(np.median(np.abs(error))),
        "fraction_with_absolute_error_below_0p2": float(np.mean(np.abs(error) < 0.2)),
    }


def load_responses(
    prefix_path: Path,
    molecule_ids: list[str],
    fraction: float,
) -> ResponseArrays:
    frame = pd.read_parquet(prefix_path)
    frame = frame.loc[
        frame.energy_group.eq("system")
        & np.isclose(frame.trajectory_fraction, fraction)
        & frame.molecule_id.astype(str).isin(molecule_ids)
    ].copy()

    def matrix(component: str, value: str) -> np.ndarray:
        selected = frame.loc[frame.component.eq(component)]
        pivoted = selected.pivot(index="molecule_id", columns="lambda", values=value)
        pivoted = pivoted.reindex(index=molecule_ids, columns=[0.1, 0.5, 0.9])
        values = pivoted.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise AssertionError(f"Missing {component}/{value} at trajectory fraction {fraction}")
        return values

    total = matrix(TOTAL_COMPONENT, "mean_kcal_mol")
    total_sem = matrix(TOTAL_COMPONENT, "five_block_sem_kcal_mol")
    component_values = np.stack(
        [matrix(component, "mean_kcal_mol") for component in COMPONENTS], axis=2
    )
    return ResponseArrays(total=total, components=component_values, total_sem=total_sem)


def fit_endpoint_ensemble(
    public_x: np.ndarray,
    public_y: np.ndarray,
    benchmark_x: np.ndarray,
    benchmark_y: np.ndarray,
    train_indices: np.ndarray,
    test_indices: np.ndarray,
) -> np.ndarray:
    x_fit = np.vstack([public_x, benchmark_x[train_indices]])
    y_fit = np.concatenate([public_y, benchmark_y[train_indices]])
    weights = np.concatenate([np.ones(len(public_y)), np.full(len(train_indices), 3.0)])
    predictions = []
    for seed in MODEL_SEEDS:
        model = endpoint_model(seed)
        model.fit(x_fit, y_fit, extratreesregressor__sample_weight=weights)
        predictions.append(model.predict(benchmark_x[test_indices]))
    return np.mean(predictions, axis=0)


def nested_endpoint_baseline(
    data,
    all_folds: np.ndarray,
    complete_indices: np.ndarray,
    cache_path: Path,
) -> np.ndarray:
    """Return [outer_fold, complete_molecule] doubly held-out SolvAI predictions."""
    if cache_path.is_file():
        return np.load(cache_path)["nested"]
    benchmark_x, public_x = data.feature_sets["F_full_solvai"]
    benchmark_y = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    public_y = data.public.delta_g_exp.to_numpy(dtype=float)
    folds = sorted(np.unique(all_folds))
    nested = np.full((len(folds), len(complete_indices)), np.nan, dtype=float)
    complete_folds = all_folds[complete_indices]
    for left, right in combinations(folds, 2):
        train = np.flatnonzero((all_folds != left) & (all_folds != right))
        for target_fold, outer_fold in ((left, right), (right, left)):
            local_test = np.flatnonzero(complete_folds == target_fold)
            global_test = complete_indices[local_test]
            nested[outer_fold, local_test] = fit_endpoint_ensemble(
                public_x,
                public_y,
                benchmark_x,
                benchmark_y,
                train,
                global_test,
            )
        print(f"  endpoint nested pair {left}/{right} complete", flush=True)
    for outer_fold in folds:
        expected = complete_folds != outer_fold
        if not np.isfinite(nested[outer_fold, expected]).all():
            raise AssertionError(f"Incomplete nested endpoint baseline for fold {outer_fold}")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, nested=nested)
    return nested


def nested_response_predictions(
    structure: np.ndarray,
    target: np.ndarray,
    folds: np.ndarray,
    cache_path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ordinary OOF and doubly held-out response predictions."""
    if cache_path.is_file():
        cache = np.load(cache_path)
        return cache["outer"], cache["nested"]
    unique_folds = sorted(np.unique(folds))
    outer = np.full_like(target, np.nan, dtype=float)
    nested = np.full((len(unique_folds), *target.shape), np.nan, dtype=float)
    for outer_fold in unique_folds:
        train = folds != outer_fold
        test = folds == outer_fold
        outer[test] = ensemble_fit_predict(structure[train], target[train], structure[test])
    for left, right in combinations(unique_folds, 2):
        train = (folds != left) & (folds != right)
        for target_fold, outer_fold in ((left, right), (right, left)):
            test = folds == target_fold
            nested[outer_fold, test] = ensemble_fit_predict(
                structure[train], target[train], structure[test]
            )
        print(f"  response nested pair {left}/{right} complete", flush=True)
    for outer_fold in unique_folds:
        expected = folds != outer_fold
        if not np.isfinite(nested[outer_fold, expected]).all():
            raise AssertionError(f"Incomplete nested response prediction for fold {outer_fold}")
    if not np.isfinite(outer).all():
        raise AssertionError("Incomplete ordinary response OOF prediction")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, outer=outer, nested=nested)
    return outer, nested


def posterior_features(
    actual: np.ndarray,
    prior_outer: np.ndarray,
    prior_nested: np.ndarray,
    sem: np.ndarray,
    folds: np.ndarray,
    outer_fold: int,
    subset: tuple[int, ...],
    conditioned: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build cross-fit outer-training and outer-test posterior curve features."""
    train_mask = folds != outer_fold
    test_mask = folds == outer_fold
    train_features = np.full((train_mask.sum(), 3), np.nan)
    train_variances = np.full((train_mask.sum(), 3), np.nan)
    train_positions = np.flatnonzero(train_mask)
    for inner_fold in sorted(np.unique(folds[train_mask])):
        inner_test = train_mask & (folds == inner_fold)
        model_rows = train_mask & (folds != inner_fold)
        if conditioned:
            model_samples = actual[model_rows] - prior_nested[outer_fold, model_rows]
            curve = fit_gaussian_curve(model_samples)
        else:
            curve = fit_gaussian_curve(actual[model_rows])
        for local_index in np.flatnonzero(inner_test):
            prior_mean = (
                prior_nested[outer_fold, local_index] + curve.mean if conditioned else curve.mean
            )
            posterior_mean, posterior_covariance = condition_gaussian(
                prior_mean,
                curve.covariance,
                subset,
                actual[local_index, subset],
                np.square(sem[local_index, subset]),
            )
            output_index = int(np.flatnonzero(train_positions == local_index)[0])
            train_features[output_index] = posterior_mean
            train_variances[output_index] = np.diag(posterior_covariance)
    if conditioned:
        curve = fit_gaussian_curve(actual[train_mask] - prior_nested[outer_fold, train_mask])
    else:
        curve = fit_gaussian_curve(actual[train_mask])
    test_features = []
    test_variances = []
    for local_index in np.flatnonzero(test_mask):
        prior_mean = prior_outer[local_index] + curve.mean if conditioned else curve.mean
        posterior_mean, posterior_covariance = condition_gaussian(
            prior_mean,
            curve.covariance,
            subset,
            actual[local_index, subset],
            np.square(sem[local_index, subset]),
        )
        test_features.append(posterior_mean)
        test_variances.append(np.diag(posterior_covariance))
    if not np.isfinite(train_features).all():
        raise AssertionError("Incomplete cross-fit posterior training features")
    return (
        train_features,
        np.asarray(test_features),
        train_variances,
        np.asarray(test_variances),
    )


def correction_prediction(
    train_features: np.ndarray,
    test_features: np.ndarray,
    baseline_train: np.ndarray,
    baseline_test: np.ndarray,
    truth_train: np.ndarray,
    train_folds: np.ndarray,
) -> tuple[np.ndarray, float]:
    residual = truth_train - baseline_train
    alpha = choose_ridge_alpha(train_features, residual, baseline_train, truth_train, train_folds)
    correction = ridge_correct(train_features, residual, test_features, alpha)
    return baseline_test + correction, alpha


def prediction_rows(
    metadata: pd.DataFrame,
    indices: np.ndarray,
    prediction: np.ndarray,
    *,
    partition: str,
    repeat: int,
    split_seed: float,
    fold: int,
    method: str,
    subset: str,
    response_scope: str,
    trajectory_fraction: float,
    alpha: float | None = None,
    shuffle_seed: int | None = None,
) -> pd.DataFrame:
    result = metadata.iloc[indices][
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            "delta_g_exp",
            "delta_g_pimd8",
        ]
    ].copy()
    result = result.rename(columns={"delta_g_exp": "y_true"})
    result["partition"] = partition
    result["repeat"] = repeat
    result["split_seed"] = split_seed
    result["fold"] = fold
    result["method"] = method
    result["lambda_subset"] = subset
    result["response_scope"] = response_scope
    result["trajectory_fraction"] = trajectory_fraction
    result["ridge_alpha"] = alpha
    result["shuffle_seed"] = shuffle_seed
    result["y_pred"] = prediction
    result["residual"] = result.y_true - result.y_pred
    result["absolute_error"] = result.residual.abs()
    return result


def paired_bootstrap(differences: np.ndarray, seed: int = BOOTSTRAP_SEED) -> dict[str, float]:
    values = np.asarray(differences, dtype=float)
    rng = np.random.default_rng(seed)
    count_less = 0
    sampled = np.empty(100_000, dtype=float)
    # Chunking avoids a large temporary integer matrix.
    offset = 0
    while offset < len(sampled):
        size = min(5000, len(sampled) - offset)
        indices = rng.integers(0, len(values), size=(size, len(values)))
        chunk = values[indices].mean(axis=1)
        sampled[offset : offset + size] = chunk
        count_less += int(np.sum(chunk < 0))
        offset += size
    low, high = np.quantile(sampled, [0.025, 0.975])
    return {
        "mean_absolute_error_difference": float(values.mean()),
        "ci_low_95": float(low),
        "ci_high_95": float(high),
        "probability_lower_error": float(count_less / len(sampled)),
        "fraction_molecules_improved": float(np.mean(values < 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, default=RELEASE_ROOT)
    parser.add_argument("--workspace-root", type=Path, default=WORKSPACE_ROOT)
    parser.add_argument("--skip-prefix-sweep", action="store_true")
    args = parser.parse_args()

    started = time.time()
    release_root = args.release_root.resolve()
    workspace_root = args.workspace_root.resolve()
    active_root = release_root / "active_solvai"
    output_root = active_root / "results/phase1"
    cache_root = output_root / "cache"
    output_root.mkdir(parents=True, exist_ok=True)
    cache_root.mkdir(parents=True, exist_ok=True)

    data = load_confirmatory_data(workspace_root, teacher_overrides(release_root))
    benchmark = data.benchmark.copy()
    benchmark_ids = benchmark.molecule_id.astype(str).tolist()
    identity = pd.read_parquet(active_root / "data/identity/probe_identity_manifest.parquet")
    identity = identity.set_index("molecule_id").reindex(benchmark_ids)
    if identity.index.isna().any() or not identity.complete_three_point_curve.notna().all():
        raise AssertionError("Probe identity manifest does not align to benchmark")
    complete_indices = np.flatnonzero(identity.complete_three_point_curve.to_numpy(bool))
    if len(complete_indices) != 72:
        raise AssertionError(f"Expected 72 complete curves, found {len(complete_indices)}")
    complete_ids = [benchmark_ids[index] for index in complete_indices]
    split_table = pd.read_csv(active_root / "data/manifests/probe_split_assignments.csv")
    split_table = split_table.loc[
        split_table.partition.isin(
            ["standardized_exclusion_primary", "standardized_exclusion_repeat"]
        )
    ].copy()
    partition_specs = []
    for (partition, repeat, split_seed), group in split_table.groupby(
        ["partition", "repeat", "split_seed"], dropna=False, sort=False
    ):
        aligned = group.set_index("molecule_id").reindex(benchmark_ids)
        partition_specs.append(
            {
                "partition": str(partition),
                "repeat": int(repeat),
                "split_seed": float(split_seed) if pd.notna(split_seed) else np.nan,
                "folds": aligned.fold.to_numpy(dtype=int),
            }
        )
    if len(partition_specs) != 6:
        raise AssertionError(f"Expected six partitions, found {len(partition_specs)}")

    stored = pd.read_parquet(
        release_root / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    prefix_path = active_root / "results/phase0/response_prefix_blocks.parquet"
    metadata = benchmark[
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            "delta_g_exp",
            "delta_g_pimd8",
        ]
    ].copy()
    truth = benchmark.delta_g_exp.to_numpy(dtype=float)
    structure = data.benchmark_structure[complete_indices]
    prediction_frames: list[pd.DataFrame] = []
    reconstruction_rows: list[dict[str, object]] = []
    integration_rows: list[pd.DataFrame] = []

    fractions = TRAJECTORY_FRACTIONS if not args.skip_prefix_sweep else (1.0,)
    responses_by_fraction = {
        fraction: load_responses(prefix_path, complete_ids, fraction) for fraction in fractions
    }

    for partition_index, spec in enumerate(partition_specs):
        partition = spec["partition"]
        repeat = spec["repeat"]
        split_seed = spec["split_seed"]
        folds_all = spec["folds"]
        folds = folds_all[complete_indices]
        label = "primary" if partition.endswith("primary") else f"repeat_{repeat}"
        print(f"\nPartition {label}", flush=True)

        stored_group = stored.loc[
            stored.partition.eq(partition)
            & stored["repeat"].eq(repeat)
            & stored.method.isin(["A_structure_only", "F_full_solvai"])
        ]
        baseline_by_method = {}
        for method in ("A_structure_only", "F_full_solvai"):
            aligned = (
                stored_group.loc[stored_group.method.eq(method)]
                .set_index("molecule_id")
                .reindex(benchmark_ids)
            )
            values = aligned.y_pred.to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise AssertionError(f"Stored parent predictions missing for {label}/{method}")
            baseline_by_method[method] = values

        nested_baseline = nested_endpoint_baseline(
            data,
            folds_all,
            complete_indices,
            cache_root / f"nested_endpoint_{label}.npz",
        )

        for outer_fold in sorted(np.unique(folds)):
            test_local = np.flatnonzero(folds == outer_fold)
            test_global = complete_indices[test_local]
            for method, output_name in (
                ("A_structure_only", "P1-0_structure_only"),
                ("F_full_solvai", "P1-A_frozen_solvai"),
            ):
                prediction_frames.append(
                    prediction_rows(
                        metadata,
                        test_global,
                        baseline_by_method[method][test_global],
                        partition=partition,
                        repeat=repeat,
                        split_seed=split_seed,
                        fold=outer_fold,
                        method=output_name,
                        subset="none",
                        response_scope="none",
                        trajectory_fraction=0.0,
                    )
                )

        for fraction in fractions:
            response = responses_by_fraction[fraction]
            target = np.column_stack([response.total, response.components.reshape(len(folds), -1)])
            predicted_outer, predicted_nested = nested_response_predictions(
                structure,
                target,
                folds,
                cache_root / f"nested_response_{label}_fraction_{fraction:.1f}.npz",
            )
            total_pred_outer = predicted_outer[:, :3]
            total_pred_nested = predicted_nested[:, :, :3]
            component_pred_outer = predicted_outer[:, 3:].reshape(len(folds), 3, 6)
            component_pred_nested = predicted_nested[:, :, 3:].reshape(5, len(folds), 3, 6)

            full_analysis = np.isclose(fraction, 1.0)
            subsets = SUBSETS if full_analysis else ((0, 1, 2),)
            scopes = ("total", "components") if full_analysis else ("total",)

            for outer_fold in sorted(np.unique(folds)):
                train_mask = folds != outer_fold
                test_mask = folds == outer_fold
                train_local = np.flatnonzero(train_mask)
                test_local = np.flatnonzero(test_mask)
                test_global = complete_indices[test_local]
                baseline_train = nested_baseline[outer_fold, train_mask]
                baseline_test = baseline_by_method["F_full_solvai"][test_global]
                truth_train = truth[complete_indices[train_local]]
                train_folds = folds[train_mask]

                for scope in scopes:
                    if scope == "total":
                        actual = response.total
                        pred_outer = total_pred_outer
                        pred_nested = total_pred_nested
                    else:
                        actual = response.components
                        pred_outer = component_pred_outer
                        pred_nested = component_pred_nested

                    for subset in subsets:
                        subset_name = key_for_subset(subset)
                        if scope == "total":
                            actual_train = actual[train_mask][:, subset]
                            actual_test = actual[test_mask][:, subset]
                            predicted_train = pred_nested[outer_fold, train_mask][:, subset]
                            predicted_test = pred_outer[test_mask][:, subset]
                        else:
                            actual_train = actual[train_mask][:, subset, :].reshape(
                                train_mask.sum(), -1
                            )
                            actual_test = actual[test_mask][:, subset, :].reshape(
                                test_mask.sum(), -1
                            )
                            predicted_train = pred_nested[outer_fold, train_mask][
                                :, subset, :
                            ].reshape(train_mask.sum(), -1)
                            predicted_test = pred_outer[test_mask][:, subset, :].reshape(
                                test_mask.sum(), -1
                            )
                        candidates = {
                            "P1-B_predicted_response": (predicted_train, predicted_test),
                            "P1-C_actual_response": (actual_train, actual_test),
                            "P1-D_actual_minus_predicted": (
                                actual_train - predicted_train,
                                actual_test - predicted_test,
                            ),
                        }
                        for method, (train_features, test_features) in candidates.items():
                            prediction, alpha = correction_prediction(
                                train_features,
                                test_features,
                                baseline_train,
                                baseline_test,
                                truth_train,
                                train_folds,
                            )
                            prediction_frames.append(
                                prediction_rows(
                                    metadata,
                                    test_global,
                                    prediction,
                                    partition=partition,
                                    repeat=repeat,
                                    split_seed=split_seed,
                                    fold=outer_fold,
                                    method=method,
                                    subset=subset_name,
                                    response_scope=scope,
                                    trajectory_fraction=fraction,
                                    alpha=alpha,
                                )
                            )

                if full_analysis:
                    # Destructive controls for the primary all-three total residual.
                    residual_train = (
                        response.total[train_mask] - total_pred_nested[outer_fold, train_mask]
                    )
                    residual_test = response.total[test_mask] - total_pred_outer[test_mask]
                    for shuffle_seed in SHUFFLE_SEEDS:
                        rng = np.random.default_rng(
                            shuffle_seed + 1000 * partition_index + outer_fold
                        )
                        shuffled_train = residual_train[rng.permutation(len(residual_train))]
                        shuffled_test = residual_test[rng.permutation(len(residual_test))]
                        prediction, alpha = correction_prediction(
                            shuffled_train,
                            shuffled_test,
                            baseline_train,
                            baseline_test,
                            truth_train,
                            train_folds,
                        )
                        prediction_frames.append(
                            prediction_rows(
                                metadata,
                                test_global,
                                prediction,
                                partition=partition,
                                repeat=repeat,
                                split_seed=split_seed,
                                fold=outer_fold,
                                method="P1-H_shuffled_residual",
                                subset="0p1_0p5_0p9",
                                response_scope="total",
                                trajectory_fraction=fraction,
                                alpha=alpha,
                                shuffle_seed=shuffle_seed,
                            )
                        )

                    # Direct numerical integration and fold-local affine calibration.
                    train_integral = integration_proxy(response.total[train_mask])
                    test_integral = integration_proxy(response.total[test_mask])
                    coefficients = fit_affine(train_integral, truth_train)
                    integrated_prediction = apply_affine(test_integral, coefficients)
                    integration_rows.append(
                        prediction_rows(
                            metadata,
                            test_global,
                            integrated_prediction,
                            partition=partition,
                            repeat=repeat,
                            split_seed=split_seed,
                            fold=outer_fold,
                            method="P1-E_actual_integral_affine",
                            subset="0p1_0p5_0p9",
                            response_scope="total",
                            trajectory_fraction=fraction,
                        ).assign(
                            integral_proxy=test_integral,
                            affine_intercept=coefficients[0],
                            affine_slope=coefficients[1],
                        )
                    )

                    for subset in SUBSETS:
                        subset_name = key_for_subset(subset)
                        for conditioned, posterior_name in (
                            (False, "P1-F_generic_posterior"),
                            (True, "P1-F_solvai_conditioned_posterior"),
                        ):
                            (
                                posterior_train,
                                posterior_test,
                                _posterior_train_variance,
                                posterior_test_variance,
                            ) = posterior_features(
                                response.total,
                                total_pred_outer,
                                total_pred_nested,
                                response.total_sem,
                                folds,
                                outer_fold,
                                subset,
                                conditioned,
                            )
                            prediction, alpha = correction_prediction(
                                posterior_train,
                                posterior_test,
                                baseline_train,
                                baseline_test,
                                truth_train,
                                train_folds,
                            )
                            prediction_frames.append(
                                prediction_rows(
                                    metadata,
                                    test_global,
                                    prediction,
                                    partition=partition,
                                    repeat=repeat,
                                    split_seed=split_seed,
                                    fold=outer_fold,
                                    method=posterior_name,
                                    subset=subset_name,
                                    response_scope="total",
                                    trajectory_fraction=fraction,
                                    alpha=alpha,
                                )
                            )

                            train_integral = integration_proxy(posterior_train)
                            test_integral = integration_proxy(posterior_test)
                            coefficients = fit_affine(train_integral, truth_train)
                            integration_rows.append(
                                prediction_rows(
                                    metadata,
                                    test_global,
                                    apply_affine(test_integral, coefficients),
                                    partition=partition,
                                    repeat=repeat,
                                    split_seed=split_seed,
                                    fold=outer_fold,
                                    method=f"{posterior_name}_integral_affine",
                                    subset=subset_name,
                                    response_scope="total",
                                    trajectory_fraction=fraction,
                                ).assign(
                                    integral_proxy=test_integral,
                                    affine_intercept=coefficients[0],
                                    affine_slope=coefficients[1],
                                )
                            )

                            hidden = sorted(set(range(3)) - set(subset))
                            for test_offset, local_index in enumerate(test_local):
                                for hidden_index in hidden:
                                    reconstruction_rows.append(
                                        {
                                            "molecule_id": complete_ids[local_index],
                                            "partition": partition,
                                            "repeat": repeat,
                                            "split_seed": split_seed,
                                            "fold": outer_fold,
                                            "posterior": posterior_name,
                                            "lambda_subset": subset_name,
                                            "hidden_lambda": float((0.1, 0.5, 0.9)[hidden_index]),
                                            "y_true": float(
                                                response.total[local_index, hidden_index]
                                            ),
                                            "y_pred": float(
                                                posterior_test[test_offset, hidden_index]
                                            ),
                                            "posterior_sd": float(
                                                np.sqrt(
                                                    max(
                                                        posterior_test_variance[
                                                            test_offset, hidden_index
                                                        ],
                                                        0.0,
                                                    )
                                                )
                                            ),
                                        }
                                    )

            print(f"  fraction {fraction:.1f} complete", flush=True)

    predictions = pd.concat(prediction_frames, ignore_index=True)
    # Mean the five shuffled predictions molecule-wise for the preregistered comparison.
    shuffled = predictions.loc[predictions.method.eq("P1-H_shuffled_residual")]
    shuffle_keys = [
        "molecule_id",
        "partition",
        "repeat",
        "split_seed",
        "fold",
        "lambda_subset",
        "response_scope",
        "trajectory_fraction",
    ]
    mean_shuffled = (
        shuffled.groupby(shuffle_keys, dropna=False, as_index=False)
        .agg(y_pred=("y_pred", "mean"))
        .merge(metadata, on="molecule_id", how="left", validate="many_to_one")
    )
    mean_shuffled = mean_shuffled.rename(columns={"delta_g_exp": "y_true"})
    mean_shuffled["method"] = "P1-H_mean_shuffled_residual"
    mean_shuffled["ridge_alpha"] = np.nan
    mean_shuffled["shuffle_seed"] = np.nan
    mean_shuffled["residual"] = mean_shuffled.y_true - mean_shuffled.y_pred
    mean_shuffled["absolute_error"] = mean_shuffled.residual.abs()
    predictions = pd.concat([predictions, mean_shuffled[predictions.columns]], ignore_index=True)
    predictions.to_parquet(output_root / "phase1_endpoint_predictions.parquet", index=False)

    metric_columns = [
        "partition",
        "repeat",
        "split_seed",
        "method",
        "lambda_subset",
        "response_scope",
        "trajectory_fraction",
        "shuffle_seed",
    ]
    metric_rows = []
    for keys, group in predictions.groupby(metric_columns, dropna=False, sort=False):
        record = dict(zip(metric_columns, keys, strict=True))
        record.update(metric_record(group.y_true.to_numpy(), group.y_pred.to_numpy()))
        record["n"] = len(group)
        metric_rows.append(record)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output_root / "phase1_endpoint_metrics.csv", index=False)

    comparisons = []
    base_keys = ["partition", "repeat", "molecule_id"]
    baseline = predictions.loc[predictions.method.eq("P1-A_frozen_solvai")].drop_duplicates(
        base_keys
    )
    baseline = baseline.set_index(base_keys)
    for keys, candidate in predictions.loc[
        ~predictions.method.isin(["P1-0_structure_only", "P1-A_frozen_solvai"])
    ].groupby(
        [
            "partition",
            "repeat",
            "split_seed",
            "method",
            "lambda_subset",
            "response_scope",
            "trajectory_fraction",
            "shuffle_seed",
        ],
        dropna=False,
        sort=False,
    ):
        aligned = candidate.set_index(base_keys)
        base = baseline.loc[aligned.index]
        difference = aligned.absolute_error.to_numpy() - base.absolute_error.to_numpy()
        comparisons.append(
            {
                "partition": keys[0],
                "repeat": keys[1],
                "split_seed": keys[2],
                "method": keys[3],
                "lambda_subset": keys[4],
                "response_scope": keys[5],
                "trajectory_fraction": keys[6],
                "shuffle_seed": keys[7],
                **paired_bootstrap(difference),
            }
        )
    comparison_frame = pd.DataFrame(comparisons)
    comparison_frame.to_csv(output_root / "phase1_paired_comparisons.csv", index=False)

    reconstruction = pd.DataFrame(reconstruction_rows)
    reconstruction["residual"] = reconstruction.y_true - reconstruction.y_pred
    reconstruction["absolute_error"] = reconstruction.residual.abs()
    for level in (0.50, 0.80, 0.90, 0.95):
        from scipy.stats import norm

        z = float(norm.ppf((1.0 + level) / 2.0))
        reconstruction[f"covered_{int(level * 100)}"] = (
            reconstruction.absolute_error <= z * reconstruction.posterior_sd
        )
        reconstruction[f"interval_width_{int(level * 100)}"] = 2.0 * z * reconstruction.posterior_sd
    reconstruction.to_parquet(
        output_root / "phase1_reconstruction_predictions.parquet", index=False
    )
    reconstruction_metrics = reconstruction.groupby(
        ["partition", "repeat", "posterior", "lambda_subset"], as_index=False
    ).agg(
        n=("absolute_error", "size"),
        mae=("absolute_error", "mean"),
        rmse=("residual", lambda value: float(np.sqrt(np.mean(np.square(value))))),
        coverage_50=("covered_50", "mean"),
        coverage_80=("covered_80", "mean"),
        coverage_90=("covered_90", "mean"),
        coverage_95=("covered_95", "mean"),
        mean_width_95=("interval_width_95", "mean"),
    )
    reconstruction_metrics.to_csv(output_root / "phase1_reconstruction_metrics.csv", index=False)

    integrations = pd.concat(integration_rows, ignore_index=True)
    integrations.to_parquet(output_root / "phase1_integration_predictions.parquet", index=False)
    integration_metric_rows = []
    for keys, group in integrations.groupby(
        ["partition", "repeat", "method", "lambda_subset"], sort=False
    ):
        record = {
            "partition": keys[0],
            "repeat": keys[1],
            "method": keys[2],
            "lambda_subset": keys[3],
            **metric_record(group.y_true.to_numpy(), group.y_pred.to_numpy()),
            "n": len(group),
        }
        pimd = group.delta_g_pimd8.to_numpy(dtype=float)
        record["mae_vs_pimd8"] = float(np.mean(np.abs(pimd - group.y_pred.to_numpy())))
        integration_metric_rows.append(record)
    pd.DataFrame(integration_metric_rows).to_csv(
        output_root / "phase1_integration_metrics.csv", index=False
    )

    # Repeated-partition primary adjudication, averaged per molecule before bootstrap.
    repeated = predictions.loc[predictions.partition.eq("standardized_exclusion_repeat")]

    def repeated_error(method: str) -> pd.Series:
        selected = repeated.loc[
            repeated.method.eq(method)
            & repeated.lambda_subset.eq("0p1_0p5_0p9")
            & repeated.response_scope.eq("total")
            & np.isclose(repeated.trajectory_fraction, 1.0)
        ]
        return selected.groupby("molecule_id").absolute_error.mean().sort_index()

    primary_error = repeated_error("P1-D_actual_minus_predicted")
    baseline_error = (
        repeated.loc[repeated.method.eq("P1-A_frozen_solvai")]
        .groupby("molecule_id")
        .absolute_error.mean()
        .sort_index()
    )
    shuffled_error = repeated_error("P1-H_mean_shuffled_residual")
    if not primary_error.index.equals(baseline_error.index) or not primary_error.index.equals(
        shuffled_error.index
    ):
        raise AssertionError("Repeated primary/bootstrap molecule sets differ")
    versus_baseline = paired_bootstrap(primary_error.to_numpy() - baseline_error.to_numpy())
    versus_shuffle = paired_bootstrap(primary_error.to_numpy() - shuffled_error.to_numpy())
    repeat_primary_rows = metrics.loc[
        metrics.partition.eq("standardized_exclusion_repeat")
        & metrics.method.eq("P1-D_actual_minus_predicted")
        & metrics.lambda_subset.eq("0p1_0p5_0p9")
        & metrics.response_scope.eq("total")
        & np.isclose(metrics.trajectory_fraction, 1.0)
        & metrics.shuffle_seed.isna()
    ].sort_values("repeat")
    improvement = -float(versus_baseline["mean_absolute_error_difference"])
    stable_sign = bool(
        (
            repeat_primary_rows.mae.to_numpy()
            < metrics.loc[
                metrics.partition.eq("standardized_exclusion_repeat")
                & metrics.method.eq("P1-A_frozen_solvai")
            ]
            .sort_values("repeat")
            .mae.to_numpy()
        ).all()
    )
    endpoint_positive = bool(
        improvement >= 0.020
        and versus_baseline["ci_high_95"] < 0
        and versus_shuffle["ci_high_95"] < 0
    )
    endpoint_negative = bool(improvement < 0.010 or not stable_sign)
    summary = {
        "status": "complete",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "freeze_commit": "a0dd986",
        "cohort_n": 72,
        "partitions": 6,
        "primary_method": "P1-D_actual_minus_predicted",
        "repeated_primary_maes": repeat_primary_rows[
            ["repeat", "mae", "rmse", "median_absolute_error"]
        ].to_dict("records"),
        "paired_vs_frozen_solvai": versus_baseline,
        "paired_vs_mean_shuffle": versus_shuffle,
        "mean_mae_improvement_vs_frozen_solvai": improvement,
        "stable_improvement_sign_across_repeats": stable_sign,
        "endpoint_gate_positive": endpoint_positive,
        "endpoint_gate_negative": endpoint_negative,
        "elapsed_seconds": time.time() - started,
        "input_hashes": {
            "freeze": sha256(active_root / "release/PHASE1_ORACLE_FREEZE.md"),
            "prefix_blocks": sha256(prefix_path),
            "splits": sha256(active_root / "data/manifests/probe_split_assignments.csv"),
            "parent_predictions": sha256(
                release_root
                / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
            ),
        },
        "output_hashes": {},
    }
    for path in sorted(output_root.glob("phase1_*")):
        if path.is_file() and path.name != "phase1_summary.json":
            summary["output_hashes"][path.name] = sha256(path)
    (output_root / "phase1_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    append_record(
        active_root / "runs/ledger.jsonl",
        {
            "run_id": "AS-P1-GATE-001",
            "stage": "phase1",
            "status": "completed",
            "command": "active_solvai/.venv/bin/python active_solvai/scripts/run_phase1_gate.py",
            "device": "CPU",
            "wall_seconds": summary["elapsed_seconds"],
            "gpu_hours": 0.0,
            "cpu_hours": summary["elapsed_seconds"] * (os.cpu_count() or 1) / 3600.0,
            "simulated_time_ps": 0.0,
            "force_evaluations": 0,
            "bead_windows": 0,
            "quality_control": "complete; see phase1_summary.json",
            "failure_reason": None,
            "freeze_commit": "a0dd986",
        },
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
