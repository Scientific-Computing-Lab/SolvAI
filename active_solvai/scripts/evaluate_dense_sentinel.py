"""Evaluate frozen replay policies on the prospective dense PIMD2 sentinels."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from active_solvai.dense import (
    condition_curve,
    curvature_order,
    curvature_scale,
    fixed_order,
    interpolate_three_point_prior,
    maximin_order,
    observed_pchip,
    rbf_covariance,
    trapezoid_weights,
    variance_reduction_order,
)

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"
CONFIG_PATH = ACTIVE_ROOT / "configs/dense_sentinel_v1.json"
LOCK_PATH = ACTIVE_ROOT / "release/DENSE_SENTINEL_CALIBRATION_LOCK.json"
CALIBRATION_PATH = ACTIVE_ROOT / "results/phase2/dense_responses_calibration.parquet"
PROSPECTIVE_PATH = ACTIVE_ROOT / "results/phase2/dense_responses_prospective.parquet"
PHASE1_PRIORS = ACTIVE_ROOT / "results/phase1/phase1_response_predictions.parquet"
OUT = ACTIVE_ROOT / "results/phase2"
INITIAL = (2, 6, 12)
Z90 = 1.6448536269514722


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bootstrap_interval(
    values: np.ndarray, seed: int, count: int, level: float
) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    sampled = values[rng.integers(0, len(values), size=(count, len(values)))].mean(axis=1)
    tail = (1.0 - level) / 2.0
    low, high = np.quantile(sampled, [tail, 1.0 - tail])
    return float(low), float(high)


def load_structure_priors(names: list[str], grid: np.ndarray) -> dict[str, np.ndarray]:
    frame = pd.read_parquet(PHASE1_PRIORS)
    frame = frame.loc[
        frame.partition.eq("standardized_exclusion_primary")
        & frame.repeat.eq(-1)
        & np.isclose(frame.trajectory_fraction, 1.0)
        & frame.component.eq("total")
        & frame.molecule_name.isin(names)
    ]
    result: dict[str, np.ndarray] = {}
    for name, group in frame.groupby("molecule_name"):
        group = group.sort_values("lambda")
        result[name] = interpolate_three_point_prior(
            group.predicted_structure_only.to_numpy(float), grid
        )
    if set(result) != set(names):
        raise AssertionError("Missing frozen structure response prior")
    return result


def timing_lookup() -> dict[tuple[str, float], float]:
    lookup: dict[tuple[str, float], float] = {}
    historical = pd.read_csv(ACTIVE_ROOT / "data/manifests/response_case_inventory.csv")
    for row in historical.loc[historical.success].to_dict("records"):
        lookup[(row["molecule_name"], float(row["lambda"]))] = float(row["elapsed_seconds"])
    status_path = ACTIVE_ROOT / "simulations/dense_pimd2/run_status_prospective.csv"
    status = pd.read_csv(status_path)
    for (name, lambda_value), group in status.groupby(["molecule_name", "lambda"]):
        lookup[(name, float(lambda_value))] = float(group.wall_seconds.sum())
    return lookup


def random_order(grid: np.ndarray, seed: int, molecule_id: str) -> list[int]:
    missing = np.array(sorted(set(range(len(grid))) - set(INITIAL)), dtype=int)
    molecule_seed = int(hashlib.sha1(molecule_id.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed + molecule_seed)
    return rng.permutation(missing).tolist()


def oracle_order(
    curve: np.ndarray,
    sem: np.ndarray,
    prior: np.ndarray,
    covariance: np.ndarray,
    weights: np.ndarray,
    noise_multiplier: float,
) -> list[int]:
    chosen = list(INITIAL)
    remaining = set(range(len(curve))) - set(chosen)
    truth = float(weights @ curve)
    order: list[int] = []
    while remaining:
        scored = []
        for candidate in remaining:
            indices = np.asarray(chosen + [candidate], dtype=int)
            posterior = condition_curve(
                prior,
                covariance,
                indices,
                curve[indices],
                noise_multiplier * sem[indices] ** 2,
            )
            estimate, _ = posterior.integral(weights)
            scored.append((abs(estimate - truth), candidate))
        _, selected = min(scored)
        chosen.append(selected)
        order.append(selected)
        remaining.remove(selected)
    return order


def replay_record(
    *,
    name: str,
    molecule_id: str,
    family: str,
    method: str,
    replicate: int,
    grid: np.ndarray,
    truth_curve: np.ndarray,
    sem: np.ndarray,
    prior: np.ndarray,
    covariance: np.ndarray,
    noise_multiplier: float,
    selected: list[int],
    elapsed_lookup: dict[tuple[str, float], float],
    direct: bool = False,
) -> dict[str, object]:
    indices = np.asarray(selected, dtype=int)
    weights = trapezoid_weights(grid)
    true_integral = -float(weights @ truth_curve)
    dense_integral_sem = float(np.sqrt(np.sum((weights * sem) ** 2)))
    if direct:
        curve = observed_pchip(grid, indices, truth_curve[indices])
        predicted_integral = -float(weights @ curve)
        integral_sd = np.nan
        covered_90 = np.nan
        width_90 = np.nan
    else:
        posterior = condition_curve(
            prior,
            covariance,
            indices,
            truth_curve[indices],
            noise_multiplier * sem[indices] ** 2,
        )
        curve = posterior.mean
        annihilation_integral, integral_sd = posterior.integral(weights)
        predicted_integral = -annihilation_integral
        covered_90 = bool(abs(predicted_integral - true_integral) <= Z90 * integral_sd)
        width_90 = 2.0 * Z90 * integral_sd
    hidden = np.array([index for index in range(len(grid)) if index not in indices], dtype=int)
    hidden_error = np.abs(curve[hidden] - truth_curve[hidden]) if len(hidden) else np.array([0.0])
    elapsed = sum(elapsed_lookup[(name, float(grid[index]))] for index in indices)
    return {
        "molecule_name": name,
        "molecule_id": molecule_id,
        "functional_group_family": family,
        "method": method,
        "schedule_replicate": replicate,
        "total_windows": len(indices),
        "observed_lambdas": ",".join(f"{grid[index]:g}" for index in indices),
        "true_integral_kcal_mol": true_integral,
        "dense_integral_sem_kcal_mol": dense_integral_sem,
        "predicted_integral_kcal_mol": predicted_integral,
        "signed_integral_error_kcal_mol": predicted_integral - true_integral,
        "absolute_integral_error_kcal_mol": abs(predicted_integral - true_integral),
        "hidden_curve_mae_kcal_mol": float(hidden_error.mean()),
        "maximum_hidden_error_kcal_mol": float(hidden_error.max()),
        "posterior_integral_sd_kcal_mol": integral_sd,
        "covered_90": covered_90,
        "interval_width_90_kcal_mol": width_90,
        "production_ps": len(indices) * 5.0,
        "bead_windows": len(indices) * 2,
        "nominal_bead_steps": len(indices) * 6250,
        "measured_window_wall_seconds": elapsed,
    }


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text())
    lock = json.loads(LOCK_PATH.read_text())
    if sha256(CALIBRATION_PATH) != lock["inputs"]["calibration_responses"]["sha256"]:
        raise AssertionError("Calibration response artifact differs from committed lock")
    grid = np.asarray(config["lambda_grid"], dtype=float)
    weights = trapezoid_weights(grid)
    prospective = pd.read_parquet(PROSPECTIVE_PATH)
    names = config["prospective_molecules"]
    counts = prospective.groupby("molecule_name").size().reindex(names)
    if not (counts == len(grid)).all():
        raise AssertionError(f"Prospective curves incomplete: {counts.to_dict()}")
    priors = load_structure_priors(names, grid)
    generic_mean = np.asarray(lock["generic_prior_mean"], dtype=float)
    generic_settings = lock["selection"]["generic"]
    solvai_settings = lock["selection"]["solvai"]
    expected_sem = np.asarray(lock["expected_sem_by_lambda"], dtype=float)
    generic_covariance = rbf_covariance(
        grid, generic_settings["amplitude"], generic_settings["lengthscale"]
    )
    elapsed_lookup = timing_lookup()
    budgets = config["evaluation_budgets_total_windows"]
    rows: list[dict[str, object]] = []
    for name in names:
        group = prospective.loc[prospective.molecule_name.eq(name)].sort_values("lambda")
        curve = group.mean_dhdl_kcal_mol.to_numpy(float)
        sem = group.five_block_sem_kcal_mol.to_numpy(float)
        molecule_id = str(group.molecule_id.iloc[0])
        family = str(group.functional_group_family.iloc[0])
        prior = priors[name]
        local = curvature_scale(prior, grid)
        solvai_covariance = rbf_covariance(
            grid,
            solvai_settings["amplitude"],
            solvai_settings["lengthscale"],
            local,
        )
        generic_expected_noise = expected_sem**2 * generic_settings["noise_inflation"]
        solvai_expected_noise = expected_sem**2 * solvai_settings["noise_inflation"]
        schedules: list[tuple[str, int, list[int], np.ndarray, np.ndarray, float, bool]] = [
            (
                "fixed_solvai_bq",
                0,
                fixed_order(grid),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                False,
            ),
            (
                "uniform_solvai_bq",
                0,
                maximin_order(grid, INITIAL),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                False,
            ),
            (
                "curvature_solvai_bq",
                0,
                curvature_order(prior, grid, INITIAL),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                False,
            ),
            (
                "generic_bq",
                0,
                variance_reduction_order(
                    generic_covariance, weights, INITIAL, generic_expected_noise
                ),
                generic_mean,
                generic_covariance,
                generic_settings["noise_inflation"],
                False,
            ),
            (
                "active_solvai_bq",
                0,
                variance_reduction_order(
                    solvai_covariance, weights, INITIAL, solvai_expected_noise
                ),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                False,
            ),
            (
                "fixed_direct",
                0,
                fixed_order(grid),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                True,
            ),
            (
                "uniform_direct",
                0,
                maximin_order(grid, INITIAL),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                True,
            ),
            (
                "oracle_non_deployable",
                0,
                oracle_order(
                    curve,
                    sem,
                    prior,
                    solvai_covariance,
                    weights,
                    solvai_settings["noise_inflation"],
                ),
                prior,
                solvai_covariance,
                solvai_settings["noise_inflation"],
                False,
            ),
        ]
        for seed in config["random_schedule_seeds"]:
            schedules.append(
                (
                    "random_solvai_bq",
                    seed,
                    random_order(grid, seed, molecule_id),
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    False,
                )
            )
        for (
            method,
            replicate,
            order,
            method_prior,
            covariance,
            noise_multiplier,
            direct,
        ) in schedules:
            for budget in budgets:
                selected = list(INITIAL) + order[: max(0, budget - len(INITIAL))]
                rows.append(
                    replay_record(
                        name=name,
                        molecule_id=molecule_id,
                        family=family,
                        method=method,
                        replicate=replicate,
                        grid=grid,
                        truth_curve=curve,
                        sem=sem,
                        prior=method_prior,
                        covariance=covariance,
                        noise_multiplier=noise_multiplier,
                        selected=selected,
                        elapsed_lookup=elapsed_lookup,
                        direct=direct,
                    )
                )
        rows.append(
            replay_record(
                name=name,
                molecule_id=molecule_id,
                family=family,
                method="full_dense_reference",
                replicate=0,
                grid=grid,
                truth_curve=curve,
                sem=sem,
                prior=prior,
                covariance=solvai_covariance,
                noise_multiplier=solvai_settings["noise_inflation"],
                selected=list(range(len(grid))),
                elapsed_lookup=elapsed_lookup,
                direct=True,
            )
        )
    predictions = pd.DataFrame(rows)
    predictions.to_parquet(OUT / "dense_replay_predictions.parquet", index=False)

    per_molecule = predictions.groupby(
        ["method", "total_windows", "molecule_id"], as_index=False
    ).agg(
        absolute_integral_error_kcal_mol=("absolute_integral_error_kcal_mol", "mean"),
        signed_integral_error_kcal_mol=("signed_integral_error_kcal_mol", "mean"),
        hidden_curve_mae_kcal_mol=("hidden_curve_mae_kcal_mol", "mean"),
        maximum_hidden_error_kcal_mol=("maximum_hidden_error_kcal_mol", "mean"),
        covered_90=("covered_90", "mean"),
        interval_width_90_kcal_mol=("interval_width_90_kcal_mol", "mean"),
        measured_window_wall_seconds=("measured_window_wall_seconds", "mean"),
    )
    metrics = per_molecule.groupby(["method", "total_windows"], as_index=False).agg(
        n=("molecule_id", "nunique"),
        integral_mae_kcal_mol=("absolute_integral_error_kcal_mol", "mean"),
        signed_integral_bias_kcal_mol=("signed_integral_error_kcal_mol", "mean"),
        hidden_curve_mae_kcal_mol=("hidden_curve_mae_kcal_mol", "mean"),
        maximum_hidden_error_kcal_mol=("maximum_hidden_error_kcal_mol", "max"),
        coverage_90=("covered_90", "mean"),
        mean_interval_width_90_kcal_mol=("interval_width_90_kcal_mol", "mean"),
        mean_measured_window_wall_seconds=("measured_window_wall_seconds", "mean"),
    )
    metrics.to_csv(OUT / "dense_replay_metrics.csv", index=False)

    comparisons: list[dict[str, object]] = []
    for budget in (5, 7):
        candidate = per_molecule.loc[
            per_molecule.method.eq("active_solvai_bq") & per_molecule.total_windows.eq(budget)
        ].set_index("molecule_id")
        for comparator in (
            "fixed_solvai_bq",
            "uniform_solvai_bq",
            "random_solvai_bq",
            "curvature_solvai_bq",
            "generic_bq",
            "fixed_direct",
            "uniform_direct",
        ):
            control = per_molecule.loc[
                per_molecule.method.eq(comparator) & per_molecule.total_windows.eq(budget)
            ].set_index("molecule_id")
            difference = (
                candidate.absolute_integral_error_kcal_mol
                - control.absolute_integral_error_kcal_mol
            ).to_numpy(float)
            low90, high90 = bootstrap_interval(
                difference,
                config["bootstrap_seed"],
                config["bootstrap_resamples"],
                0.90,
            )
            low95, high95 = bootstrap_interval(
                difference,
                config["bootstrap_seed"],
                config["bootstrap_resamples"],
                0.95,
            )
            comparisons.append(
                {
                    "total_windows": budget,
                    "candidate": "active_solvai_bq",
                    "comparator": comparator,
                    "candidate_minus_comparator_mae": float(difference.mean()),
                    "ci90_low": low90,
                    "ci90_high": high90,
                    "ci95_low": low95,
                    "ci95_high": high95,
                    "fraction_candidate_improved": float(np.mean(difference < 0)),
                }
            )
    comparison_frame = pd.DataFrame(comparisons)
    comparison_frame.to_csv(OUT / "dense_replay_paired_comparisons.csv", index=False)
    print(metrics.to_string(index=False))
    print(comparison_frame.to_string(index=False))


if __name__ == "__main__":
    main()
