#!/usr/bin/env python3
"""Run the prospectively frozen zero-new-simulation identifiability gate."""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import UndefinedMetricWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from active_solvai_v3.diagnostics import complementary_log_difficulty, prefix_diagnostics

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "active_solvai_v3/data/derived"
OUT = ROOT / "active_solvai_v3/results/gate1"
PREFIX_PS = (0.5, 1.0, 2.0, 3.0)
STABILIZERS = (0.1, 0.25, 0.5)
ALPHAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
SHUFFLE_SEEDS = (88031, 88032, 88033, 88034, 88035)
BOOTSTRAP_SEED = 20260903
BOOTSTRAP_REPLICATES = 100_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def trapezoid_weights(grid: np.ndarray) -> np.ndarray:
    grid = np.asarray(grid, dtype=float)
    if grid.ndim != 1 or len(grid) < 2 or np.any(np.diff(grid) <= 0):
        raise ValueError("lambda grid must be strictly increasing")
    weights = np.zeros_like(grid)
    differences = np.diff(grid)
    weights[:-1] += differences / 2
    weights[1:] += differences / 2
    return weights


def build_diagnostic_table(time_series: pd.DataFrame, molecule_features: pd.DataFrame) -> pd.DataFrame:
    grid = np.sort(time_series["lambda"].unique())
    weight_map = dict(zip(grid, trapezoid_weights(grid), strict=True))
    rows: list[dict[str, object]] = []
    for (molecule_id, lambda_value), trajectory in time_series.groupby(
        ["molecule_id", "lambda"], sort=True
    ):
        trajectory = trajectory.sort_values("frame_index")
        if trajectory.frame_index.tolist() != list(range(1, 51)):
            raise AssertionError(f"non-canonical frames for {molecule_id}/{lambda_value}")
        block_four = trajectory.loc[trajectory.block_1ps == 4, "dhdl_kcal_mol"].to_numpy()
        block_five = trajectory.loc[trajectory.block_1ps == 5, "dhdl_kcal_mol"].to_numpy()
        targets = {
            stabilizer: complementary_log_difficulty(block_four, block_five, stabilizer)
            for stabilizer in STABILIZERS
        }
        for prefix_ps in PREFIX_PS:
            count = int(round(prefix_ps / 0.1))
            prefix = trajectory.iloc[:count]
            diagnostics = prefix_diagnostics(
                prefix.dhdl_kcal_mol.to_numpy(), prefix.time_ps.to_numpy()
            )
            row: dict[str, object] = {
                "molecule_id": molecule_id,
                "molecule_name": trajectory.molecule_name.iloc[0],
                "canonical_smiles": trajectory.canonical_smiles.iloc[0],
                "lambda": float(lambda_value),
                "lambda_squared": float(lambda_value**2),
                "lambda_centered": float(lambda_value - 0.5),
                "lambda_centered_squared": float((lambda_value - 0.5) ** 2),
                "endpoint_indicator": int(lambda_value in (0.0, 1.0)),
                "trapezoid_weight": float(weight_map[lambda_value]),
                "prefix_ps": prefix_ps,
                "prefix_temperature_mean": float(prefix.temperature_k.mean()),
                "prefix_density_mean": float(prefix.density_g_cm3.mean()),
                **diagnostics,
            }
            for stabilizer, target in targets.items():
                row[f"target_log_difficulty_s{str(stabilizer).replace('.', 'p')}"] = target
            rows.append(row)
    result = pd.DataFrame(rows).merge(molecule_features, on=["molecule_id", "molecule_name", "canonical_smiles"], how="left")
    if result.isna().all(axis=0).any():
        raise AssertionError("Gate-1 table contains an entirely missing column")
    return add_interactions(result)


def add_interactions(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    molecule_columns = [
        column
        for column in result
        if column.startswith("structure_") or column.startswith("response__")
    ]
    for column in molecule_columns:
        result[f"{column}__x_lambda_centered"] = result[column] * result["lambda_centered"]
        result[f"{column}__x_lambda_centered_squared"] = (
            result[column] * result["lambda_centered_squared"]
        )
    return result


def feature_sets(frame: pd.DataFrame) -> dict[str, list[str]]:
    protocol = ["lambda", "lambda_squared", "endpoint_indicator", "trapezoid_weight"]
    structure = [column for column in frame if column.startswith("structure_")]
    response = [column for column in frame if column.startswith("response__")]
    observed = [
        "prefix_mean",
        "prefix_variance",
        "naive_sem",
        "batch_variance_rate_2",
        "batch_variance_rate_5",
        "overlap_batch_variance_rate",
        "newey_west_variance_rate",
        "lag1_autocorrelation",
        "iat_initial_positive",
        "effective_sample_size",
        "half_mean_difference",
        "half_variance_ratio_log",
        "half_ks_distance",
        "linear_drift_per_ps",
        "rank_drift",
        "unresolved_equilibration",
        "prefix_temperature_mean",
        "prefix_density_mean",
    ]
    return {
        "lambda_protocol": protocol,
        "structure_cold_start": [*protocol, *structure],
        "solvai_cold_start": [*protocol, *response],
        "generic_observed": [*protocol, *observed],
        "structure_conditioned": [*protocol, *observed, *structure],
        "solvai_conditioned": [*protocol, *observed, *response],
    }


def pipeline(alpha: float):
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
        StandardScaler(),
        Ridge(alpha=alpha),
    )


def select_alpha(frame: pd.DataFrame, features: list[str], target: str) -> float:
    groups = frame.molecule_id.to_numpy()
    splitter = LeaveOneGroupOut()
    scores: dict[float, float] = {}
    for alpha in ALPHAS:
        predictions = np.empty(len(frame), dtype=float)
        for train_index, test_index in splitter.split(frame, groups=groups):
            model = pipeline(alpha)
            model.fit(frame.iloc[train_index][features], frame.iloc[train_index][target])
            predictions[test_index] = model.predict(frame.iloc[test_index][features])
        scores[alpha] = float(mean_absolute_error(frame[target], predictions))
    minimum = min(scores.values())
    return max(alpha for alpha, score in scores.items() if np.isclose(score, minimum))


def shuffle_outer_training_responses(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    result = frame.copy()
    response_main = [
        column
        for column in result
        if column.startswith("response__") and "__x_lambda" not in column
    ]
    molecules = sorted(result.molecule_id.unique())
    shuffled = np.random.default_rng(seed).permutation(molecules)
    mapping = dict(zip(molecules, shuffled, strict=True))
    lookup = result[["molecule_id", *response_main]].drop_duplicates("molecule_id").set_index("molecule_id")
    for molecule, source in mapping.items():
        mask = result.molecule_id == molecule
        result.loc[mask, response_main] = lookup.loc[source, response_main].to_numpy()
    for column in response_main:
        result[f"{column}__x_lambda_centered"] = result[column] * result["lambda_centered"]
        result[f"{column}__x_lambda_centered_squared"] = (
            result[column] * result["lambda_centered_squared"]
        )
    return result


def cross_validated_predictions(
    frame: pd.DataFrame,
    features: list[str],
    target: str,
    *,
    shuffle_seed: int | None = None,
) -> pd.DataFrame:
    outputs: list[pd.DataFrame] = []
    groups = frame.molecule_id.to_numpy()
    for outer_train, outer_test in LeaveOneGroupOut().split(frame, groups=groups):
        train = frame.iloc[outer_train].copy()
        test = frame.iloc[outer_test].copy()
        if shuffle_seed is not None:
            train = shuffle_outer_training_responses(train, shuffle_seed)
        alpha = select_alpha(train, features, target)
        model = pipeline(alpha)
        model.fit(train[features], train[target])
        output = test[["molecule_id", "molecule_name", "lambda", "prefix_ps"]].copy()
        output["target"] = test[target].to_numpy(float)
        output["prediction"] = model.predict(test[features])
        output["selected_alpha"] = alpha
        outputs.append(output)
    return pd.concat(outputs, ignore_index=True)


def model_metrics(predictions: pd.DataFrame) -> dict[str, float]:
    errors = np.abs(predictions.prediction - predictions.target)
    correlations: list[float] = []
    top_labels: list[np.ndarray] = []
    top_scores: list[np.ndarray] = []
    for _, part in predictions.groupby("molecule_id"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            correlation = spearmanr(part.target, part.prediction).statistic
        if np.isfinite(correlation):
            correlations.append(float(correlation))
        threshold = float(part.target.quantile(0.75))
        top_labels.append((part.target >= threshold).astype(int).to_numpy())
        top_scores.append(part.prediction.to_numpy(float))
    labels = np.concatenate(top_labels)
    scores = np.concatenate(top_scores)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UndefinedMetricWarning)
        auc = roc_auc_score(labels, scores)
    return {
        "mae": float(errors.mean()),
        "rmse": float(np.sqrt(np.mean((predictions.prediction - predictions.target) ** 2))),
        "mean_within_molecule_spearman": float(np.mean(correlations)),
        "top_quartile_auroc": float(auc),
    }


def paired_bootstrap(differences: pd.Series) -> dict[str, float]:
    values = differences.to_numpy(float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(0, len(values), size=(BOOTSTRAP_REPLICATES, len(values)))
    means = values[draws].mean(axis=1)
    return {
        "mean_difference": float(values.mean()),
        "ci90_low": float(np.quantile(means, 0.05)),
        "ci90_high": float(np.quantile(means, 0.95)),
        "ci95_low": float(np.quantile(means, 0.025)),
        "ci95_high": float(np.quantile(means, 0.975)),
        "fraction_molecules_below_zero": float(np.mean(values < 0)),
    }


def molecule_errors(predictions: pd.DataFrame) -> pd.Series:
    return predictions.assign(error=lambda value: abs(value.prediction - value.target)).groupby(
        "molecule_id"
    ).error.mean()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    time_path = DATA / "gate1_time_series.parquet"
    feature_path = DATA / "gate1_molecule_features.parquet"
    time_series = pd.read_parquet(time_path)
    molecules = pd.read_parquet(feature_path)
    diagnostics = build_diagnostic_table(time_series, molecules)
    diagnostic_path = OUT / "gate1_prefix_diagnostics.parquet"
    diagnostics.to_parquet(diagnostic_path, index=False)
    features_by_model = feature_sets(diagnostics)

    prediction_frames: list[pd.DataFrame] = []
    for stabilizer in STABILIZERS:
        target = f"target_log_difficulty_s{str(stabilizer).replace('.', 'p')}"
        for prefix_ps in PREFIX_PS:
            subset = diagnostics[np.isclose(diagnostics.prefix_ps, prefix_ps)].reset_index(drop=True)
            for model_name, columns in features_by_model.items():
                result = cross_validated_predictions(subset, columns, target)
                result["model"] = model_name
                result["stabilizer"] = stabilizer
                prediction_frames.append(result)
            for seed in SHUFFLE_SEEDS:
                result = cross_validated_predictions(
                    subset,
                    features_by_model["solvai_conditioned"],
                    target,
                    shuffle_seed=seed,
                )
                result["model"] = f"solvai_conditioned_shuffled_{seed}"
                result["stabilizer"] = stabilizer
                prediction_frames.append(result)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    prediction_path = OUT / "gate1_oof_predictions.parquet"
    predictions.to_parquet(prediction_path, index=False)

    metric_rows: list[dict[str, object]] = []
    for (prefix_ps, stabilizer, model), part in predictions.groupby(
        ["prefix_ps", "stabilizer", "model"], sort=True
    ):
        metric_rows.append(
            {
                "prefix_ps": prefix_ps,
                "stabilizer": stabilizer,
                "model": model,
                **model_metrics(part),
            }
        )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "gate1_model_metrics.csv", index=False)

    comparison_pairs = [
        ("solvai_conditioned", "generic_observed", "solvai_minus_generic"),
        ("generic_observed", "lambda_protocol", "generic_minus_lambda"),
        ("structure_conditioned", "generic_observed", "structure_minus_generic"),
        ("solvai_conditioned", "structure_conditioned", "solvai_minus_structure"),
        ("solvai_cold_start", "lambda_protocol", "solvai_cold_minus_lambda"),
        ("solvai_cold_start", "structure_cold_start", "solvai_cold_minus_structure"),
    ]
    comparison_rows: list[dict[str, object]] = []
    for (prefix_ps, stabilizer), part in predictions.groupby(["prefix_ps", "stabilizer"]):
        model_parts = {name: group for name, group in part.groupby("model")}
        for left, right, name in comparison_pairs:
            left_errors = molecule_errors(model_parts[left])
            right_errors = molecule_errors(model_parts[right])
            comparison_rows.append(
                {
                    "prefix_ps": prefix_ps,
                    "stabilizer": stabilizer,
                    "comparison": name,
                    **paired_bootstrap(left_errors - right_errors),
                }
            )
        aligned = molecule_errors(model_parts["solvai_conditioned"])
        shuffled_errors = pd.concat(
            [
                molecule_errors(model_parts[f"solvai_conditioned_shuffled_{seed}"]).rename(seed)
                for seed in SHUFFLE_SEEDS
            ],
            axis=1,
        ).mean(axis=1)
        comparison_rows.append(
            {
                "prefix_ps": prefix_ps,
                "stabilizer": stabilizer,
                "comparison": "aligned_minus_mean_shuffled",
                **paired_bootstrap(aligned - shuffled_errors),
            }
        )
    comparisons = pd.DataFrame(comparison_rows)
    comparisons.to_csv(OUT / "gate1_paired_comparisons.csv", index=False)

    reliability_rows: list[dict[str, object]] = []
    diagnostic_columns = feature_sets(diagnostics)["generic_observed"][4:]
    for prefix_ps in PREFIX_PS:
        part = diagnostics[np.isclose(diagnostics.prefix_ps, prefix_ps)]
        target = part["target_log_difficulty_s0p25"]
        for column in diagnostic_columns:
            values = pd.to_numeric(part[column], errors="coerce")
            finite = np.isfinite(values) & np.isfinite(target)
            correlation = spearmanr(values[finite], target[finite]).statistic if finite.sum() >= 3 else np.nan
            reliability_rows.append(
                {
                    "prefix_ps": prefix_ps,
                    "diagnostic": column,
                    "finite_fraction": float(finite.mean()),
                    "spearman_with_complementary_log_difficulty": float(correlation),
                }
            )
    reliability = pd.DataFrame(reliability_rows)
    reliability.to_csv(OUT / "gate1_diagnostic_reliability.csv", index=False)

    primary_metrics = metrics[
        np.isclose(metrics.prefix_ps, 1.0) & np.isclose(metrics.stabilizer, 0.25)
    ].set_index("model")
    primary_comparisons = comparisons[
        np.isclose(comparisons.prefix_ps, 1.0) & np.isclose(comparisons.stabilizer, 0.25)
    ].set_index("comparison")
    two_ps = comparisons[
        np.isclose(comparisons.prefix_ps, 2.0) & np.isclose(comparisons.stabilizer, 0.25)
    ].set_index("comparison")
    solvai_mae = float(primary_metrics.loc["solvai_conditioned", "mae"])
    generic_mae = float(primary_metrics.loc["generic_observed", "mae"])
    shuffled_maes = [
        float(primary_metrics.loc[f"solvai_conditioned_shuffled_{seed}", "mae"])
        for seed in SHUFFLE_SEEDS
    ]
    solvai_relative_gain = (generic_mae - solvai_mae) / generic_mae
    shuffle_relative_gain = (float(np.mean(shuffled_maes)) - solvai_mae) / float(
        np.mean(shuffled_maes)
    )
    solvai_comparison = primary_comparisons.loc["solvai_minus_generic"]
    shuffle_comparison = primary_comparisons.loc["aligned_minus_mean_shuffled"]
    ai_gate_passed = bool(
        solvai_relative_gain >= 0.10
        and solvai_comparison.ci90_high < 0
        and shuffle_relative_gain >= 0.10
        and shuffle_comparison.ci90_high < 0
        and primary_metrics.loc["solvai_conditioned", "mean_within_molecule_spearman"] >= 0.30
        and two_ps.loc["solvai_minus_generic", "mean_difference"] < 0
    )
    generic_mae = float(primary_metrics.loc["generic_observed", "mae"])
    lambda_mae = float(primary_metrics.loc["lambda_protocol", "mae"])
    generic_comparison = primary_comparisons.loc["generic_minus_lambda"]
    generic_gate_passed = bool(
        (lambda_mae - generic_mae) / lambda_mae >= 0.10 and generic_comparison.ci90_high < 0
    )
    canonical = {
        "schema_version": 1,
        "protocol_commit": "4608c21",
        "n_molecules": int(diagnostics.molecule_id.nunique()),
        "n_windows": int(
            diagnostics[np.isclose(diagnostics.prefix_ps, 1.0)].shape[0]
        ),
        "primary_prefix_ps": 1.0,
        "primary_stabilizer_kcal_mol": 0.25,
        "generic_gate_passed": generic_gate_passed,
        "solvai_identifiability_gate_passed": ai_gate_passed,
        "primary_metrics": primary_metrics.reset_index().to_dict(orient="records"),
        "primary_comparisons": primary_comparisons.reset_index().to_dict(orient="records"),
        "solvai_relative_mae_gain_vs_generic": solvai_relative_gain,
        "solvai_relative_mae_gain_vs_mean_shuffled": shuffle_relative_gain,
        "independence_limit": (
            "Complementary blocks of one 5-ps trajectory; not independent-replica validation."
        ),
    }
    canonical_path = OUT / "gate1_canonical_metrics.json"
    canonical_path.write_text(json.dumps(canonical, indent=2) + "\n")

    outputs = [
        diagnostic_path,
        prediction_path,
        OUT / "gate1_model_metrics.csv",
        OUT / "gate1_paired_comparisons.csv",
        OUT / "gate1_diagnostic_reliability.csv",
        canonical_path,
    ]
    manifest = {
        "schema_version": 1,
        "inputs": {str(path): sha256(path) for path in (time_path, feature_path)},
        "outputs": {str(path): sha256(path) for path in outputs},
    }
    (OUT / "gate1_artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(canonical, indent=2))


if __name__ == "__main__":
    main()
