"""Strict structure-only distillation of the ARROW/PIMD response hierarchy.

The response student learns ``(molecule, lambda) -> dH/dlambda`` from short
PIMD2 teachers.  Predicted curves are integrated and passed to the final
experimental head.  In every outer fold, measured physics for outer-test
molecules is excluded; inner-OOF physics predictions are used for outer-train
molecules.  Consequently every final prediction is structure-only.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from arrow_distill.experiments import REGIMES, mae, make_prediction_rows
from run_nested_smd_teacher_confirmation import fit_model, load_problem
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline

LAMBDA_VALUES = np.asarray([0.1, 0.5, 0.9], dtype=np.float32)
RESPONSE_COMPONENTS = (
    "lig_slv__dhdl_mean",
    "lig_slv__dhdl_coul_mean",
    "lig_slv__dhdl_vdw_mean",
    "lig__dhdl_pol_mean",
)


def teacher_table(benchmark: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    table = pd.read_parquet(ROOT / "results/pimd2_multilambda_teacher.parquet")
    aligned = (
        table.drop_duplicates("molecule_id")
        .set_index("molecule_id")
        .reindex(benchmark.molecule_id.astype(str))
    )
    response = np.full(
        (len(benchmark), len(LAMBDA_VALUES), len(RESPONSE_COMPONENTS)),
        np.nan,
        dtype=np.float32,
    )
    for lambda_index, code in enumerate(("01", "05", "09")):
        success = aligned[f"success__lambda{code}"].fillna(False).to_numpy(bool)
        for component_index, component in enumerate(RESPONSE_COMPONENTS):
            values = pd.to_numeric(aligned[f"{component}__lambda{code}"], errors="coerce").to_numpy(
                float, copy=True
            )
            values[~success] = np.nan
            response[:, lambda_index, component_index] = values
    hierarchy = benchmark[["delta_g_classical_arrow", "nqe_residual", "delta_g_pimd8"]].to_numpy(
        dtype=np.float32
    )
    return response, hierarchy


def lambda_design(x: np.ndarray, lambda_values: np.ndarray) -> np.ndarray:
    molecule = np.repeat(x, len(lambda_values), axis=0)
    values = np.tile(lambda_values, len(x)).reshape(-1, 1)
    basis = np.column_stack([values, values**2, values**3])
    return np.column_stack([molecule, basis]).astype(np.float32)


def response_model(seed: int, trees: int) -> object:
    return make_pipeline(
        SimpleImputer(strategy="median"),
        ExtraTreesRegressor(
            n_estimators=trees,
            min_samples_leaf=3,
            max_features=0.7,
            n_jobs=8,
            random_state=seed,
        ),
    )


def fit_response(
    x: np.ndarray,
    response: np.ndarray,
    molecule_indices: np.ndarray,
    *,
    seed: int,
    trees: int,
) -> object:
    x_rows = lambda_design(x[molecule_indices], LAMBDA_VALUES)
    y_rows = response[molecule_indices].reshape(-1, len(RESPONSE_COMPONENTS))
    valid = np.isfinite(y_rows).all(axis=1)
    estimator = response_model(seed, trees)
    estimator.fit(x_rows[valid], y_rows[valid])
    return estimator


def predict_response(estimator: object, x: np.ndarray) -> np.ndarray:
    prediction = estimator.predict(lambda_design(x, LAMBDA_VALUES))
    return prediction.reshape(len(x), len(LAMBDA_VALUES), len(RESPONSE_COMPONENTS))


def curve_latent(response: np.ndarray) -> np.ndarray:
    y01, y05, y09 = response[:, 0], response[:, 1], response[:, 2]
    integral = 0.1 * y01 + (0.8 / 6.0) * (y01 + 4.0 * y05 + y09) + 0.1 * y09
    slope = (y09 - y01) / 0.8
    curvature = y01 - 2.0 * y05 + y09
    return np.column_stack([integral, slope, curvature]).astype(np.float32)


def hierarchy_model(seed: int, trees: int) -> object:
    return make_pipeline(
        SimpleImputer(strategy="median"),
        ExtraTreesRegressor(
            n_estimators=trees,
            min_samples_leaf=3,
            max_features=0.7,
            n_jobs=8,
            random_state=seed,
        ),
    )


def fit_hierarchy(
    x: np.ndarray,
    hierarchy: np.ndarray,
    indices: np.ndarray,
    *,
    seed: int,
    trees: int,
) -> object:
    valid = np.isfinite(hierarchy[indices]).all(axis=1)
    estimator = hierarchy_model(seed, trees)
    estimator.fit(x[indices][valid], hierarchy[indices][valid])
    return estimator


def distill_outer_fold(
    x_benchmark: np.ndarray,
    x_public: np.ndarray,
    response: np.ndarray,
    hierarchy: np.ndarray,
    outer_train: np.ndarray,
    outer_test: np.ndarray,
    *,
    seed: int,
    trees: int,
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], np.ndarray]:
    response_train = np.full(
        (len(outer_train), len(LAMBDA_VALUES), len(RESPONSE_COMPONENTS)), np.nan
    )
    hierarchy_train = np.full((len(outer_train), hierarchy.shape[1]), np.nan)
    inner = KFold(n_splits=5, shuffle=True, random_state=seed)
    for inner_fold, (fit_position, held_position) in enumerate(inner.split(outer_train)):
        fit_indices = outer_train[fit_position]
        held_indices = outer_train[held_position]
        response_estimator = fit_response(
            x_benchmark,
            response,
            fit_indices,
            seed=seed + inner_fold,
            trees=trees,
        )
        response_train[held_position] = predict_response(
            response_estimator, x_benchmark[held_indices]
        )
        hierarchy_estimator = fit_hierarchy(
            x_benchmark,
            hierarchy,
            fit_indices,
            seed=seed + 101 + inner_fold,
            trees=trees,
        )
        hierarchy_train[held_position] = hierarchy_estimator.predict(x_benchmark[held_indices])

    response_estimator = fit_response(
        x_benchmark, response, outer_train, seed=seed + 1009, trees=trees
    )
    hierarchy_estimator = fit_hierarchy(
        x_benchmark, hierarchy, outer_train, seed=seed + 2017, trees=trees
    )
    response_test = predict_response(response_estimator, x_benchmark[outer_test])
    response_public = predict_response(response_estimator, x_public)
    hierarchy_test = hierarchy_estimator.predict(x_benchmark[outer_test])
    hierarchy_public = hierarchy_estimator.predict(x_public)
    blocks = {
        "response": (
            curve_latent(response_train),
            curve_latent(response_test),
            curve_latent(response_public),
        ),
        "hierarchy": (hierarchy_train, hierarchy_test, hierarchy_public),
    }
    return blocks, response_test


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regimes", nargs="+", choices=tuple(REGIMES), default=["random_oof"])
    parser.add_argument("--physics-trees", type=int, default=180)
    parser.add_argument("--head-trees", type=int, default=360)
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 29, 47])
    args = parser.parse_args()

    benchmark, public, feature_sets, source_mask = load_problem(
        no_static_arrow=True,
        confsolv_response=True,
    )
    x_benchmark, x_public = feature_sets["narrow response + SMD + ConfSolv response"]
    truth = benchmark.delta_g_exp.to_numpy(float)
    public_truth = public.delta_g_exp.to_numpy(float)[source_mask]
    response, hierarchy = teacher_table(benchmark)
    if np.isfinite(response).all(axis=(1, 2)).sum() < 50:
        raise SystemExit("Fewer than 50 complete multi-lambda PIMD2 teachers")

    groups = {
        "A: structure/response baseline": (),
        "B1: +distilled classical-NQE-PIMD hierarchy": ("hierarchy",),
        "B2: +distilled PIMD2 lambda response": ("response",),
        "B: +full distilled physics hierarchy": ("hierarchy", "response"),
    }
    integrated_method = "C: integrated predicted dH/dlambda (affine calibration)"
    rows: list[pd.DataFrame] = []
    curve_rows: list[dict[str, object]] = []
    for regime in args.regimes:
        folds = benchmark[REGIMES[regime]].to_numpy()
        predictions = {
            name: np.full(len(benchmark), np.nan) for name in (*groups, integrated_method)
        }
        for fold in sorted(np.unique(folds)):
            outer_test = np.flatnonzero(folds == fold)
            outer_train = np.flatnonzero(folds != fold)
            blocks, response_test = distill_outer_fold(
                x_benchmark,
                x_public,
                response,
                hierarchy,
                outer_train,
                outer_test,
                seed=20260826 + int(fold),
                trees=args.physics_trees,
            )
            # The first response-latent column is the explicitly integrated
            # total dH/dlambda curve.  Calibration sees only inner-OOF curve
            # predictions and experimental labels from the outer-training set.
            integrated_train = blocks["response"][0][:, [0]]
            integrated_test = blocks["response"][1][:, [0]]
            integration_calibration = LinearRegression().fit(integrated_train, truth[outer_train])
            predictions[integrated_method][outer_test] = integration_calibration.predict(
                integrated_test
            )
            for local_index, molecule_index in enumerate(outer_test):
                for lambda_index, lambda_value in enumerate(LAMBDA_VALUES):
                    for component_index, component in enumerate(RESPONSE_COMPONENTS):
                        curve_rows.append(
                            {
                                "molecule_id": benchmark.iloc[molecule_index].molecule_id,
                                "molecule_name": benchmark.iloc[molecule_index].molecule_name,
                                "regime": regime,
                                "fold": fold,
                                "lambda": float(lambda_value),
                                "component": component,
                                "y_true": float(
                                    response[molecule_index, lambda_index, component_index]
                                ),
                                "y_pred": float(
                                    response_test[local_index, lambda_index, component_index]
                                ),
                            }
                        )
            for group_name, selected_blocks in groups.items():
                fold_predictions = []
                train_parts = [x_benchmark[outer_train]]
                test_parts = [x_benchmark[outer_test]]
                public_parts = [x_public]
                for selected in selected_blocks:
                    train_block, test_block, public_block = blocks[selected]
                    train_parts.append(train_block)
                    test_parts.append(test_block)
                    public_parts.append(public_block)
                train_x = np.column_stack(train_parts)
                test_x = np.column_stack(test_parts)
                public_x = np.column_stack(public_parts)
                for model_seed in args.seeds:
                    estimator = fit_model(
                        public_x[source_mask],
                        public_truth,
                        train_x,
                        truth[outer_train],
                        seed=model_seed,
                        trees=args.head_trees,
                    )
                    fold_predictions.append(estimator.predict(test_x))
                predictions[group_name][outer_test] = np.mean(fold_predictions, axis=0)
            print(regime, fold, flush=True)
        for name, prediction in predictions.items():
            rows.append(
                make_prediction_rows(
                    benchmark,
                    prediction,
                    f"Multi-lambda physics distillation {name}",
                    regime,
                    False,
                    details=["outer-test physics excluded; structure-only inference"]
                    * len(benchmark),
                )
            )
            print(regime, name, mae(truth, prediction), flush=True)

    output = pd.concat(rows, ignore_index=True)
    output_path = ROOT / "results/multilambda_physics_distillation_oof.parquet"
    output.to_parquet(output_path, index=False)
    metrics = (
        output.groupby(["method", "regime"], as_index=False)
        .agg(n=("absolute_error", "size"), mae=("absolute_error", "mean"))
        .sort_values(["regime", "mae"])
    )
    metrics.to_csv(output_path.with_suffix(".csv"), index=False)
    curve = pd.DataFrame(curve_rows)
    curve["absolute_error"] = (curve.y_pred - curve.y_true).abs()
    curve.to_parquet(ROOT / "results/multilambda_response_student_oof.parquet", index=False)
    curve_metrics_rows = []
    for keys, frame in curve[np.isfinite(curve.y_true)].groupby(
        ["regime", "lambda", "component"], sort=False
    ):
        target_std = float(frame.y_true.std(ddof=0))
        curve_metrics_rows.append(
            {
                "regime": keys[0],
                "lambda": keys[1],
                "component": keys[2],
                "n": len(frame),
                "mae": float(frame.absolute_error.mean()),
                "rmse": float(np.sqrt(np.mean(np.square(frame.y_pred - frame.y_true)))),
                "pearson_r": float(frame.y_true.corr(frame.y_pred)),
                "target_std": target_std,
                "normalized_mae": (
                    float(frame.absolute_error.mean() / target_std) if target_std > 0 else np.nan
                ),
            }
        )
    curve_metrics = pd.DataFrame(curve_metrics_rows)
    curve_metrics.to_csv(ROOT / "results/multilambda_response_student_metrics.csv", index=False)
    metadata = {
        "inference": "SMILES/RDKit structure only; no probe, trajectory, or static ARROW input",
        "teacher_use": "outer-training folds only; inner-OOF predictions for experimental head",
        "lambda_values": LAMBDA_VALUES.tolist(),
        "response_components": list(RESPONSE_COMPONENTS),
        "complete_teacher_molecules": int(np.isfinite(response).all(axis=(1, 2)).sum()),
    }
    output_path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(metrics.to_string(index=False))
    print(curve_metrics.to_string(index=False))


if __name__ == "__main__":
    main()
