"""Run the frozen cross-block oracle-headroom diagnostic.

The protocol is committed in ORACLE_INDEPENDENT_NOISE_DIAGNOSTIC_FREEZE.md.
This script performs no fitting and launches no simulation.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from active_solvai.dense import (
    condition_curve,
    curvature_scale,
    fixed_order,
    interpolate_three_point_prior,
    maximin_order,
    observed_pchip,
    rbf_covariance,
    trapezoid_weights,
    variance_reduction_order,
)
from active_solvai.ledger import append_record
from active_solvai.probes import read_energy

ROOT = Path(__file__).resolve().parents[2]
ACTIVE = ROOT / "active_solvai"
PHASE2 = ACTIVE / "results/phase2"
OUT = ACTIVE / "results/v2_diagnostics/oracle_independent_noise"
FIGURES = ACTIVE / "figures/v2_diagnostics"
REPORT = ACTIVE / "reports/ORACLE_INDEPENDENT_NOISE_DIAGNOSTIC.md"
INITIAL = (2, 6, 12)
Z90 = 1.6448536269514722
BOOTSTRAP_SEED = 20260903
BOOTSTRAP_RESAMPLES = 100_000
EXPECTED_HASHES = {
    "responses": "4ce191a2ebfabb10dfcf5e11f98fef91b0bd5b2993f066ab47cfa9cc261c1097",
    "lock": "3d1f411c5320bf17b0bf83d7f6280645848136c88fed7d7be1a54e9e0fd4975a",
    "config": "fc6afd2a4f8a255a641228df624769e6ed437dc4cc69dc244480a33c2936594d",
    "priors": "9d9a514e2cabec3ea6f43600d9987d886ca11f2c96e2cc4f5f0b18f79397c91b",
    "manifest": "9c46ad5b78af027f3c82d84ff300c1ea839fbd13650eba704bbfef95f5871358",
    "energy_records": "e0585cc8e1f0d19754b43eec6eaa54c760fcbf106fe76315e46fbf3b7374451e",
}
SPLITS = {
    0: ((0, 1), (2, 3)),
    1: ((0, 2), (1, 3)),
    2: ((0, 3), (1, 2)),
    3: ((1, 2), (0, 3)),
    4: ((1, 3), (0, 2)),
    5: ((2, 3), (0, 1)),
}
REVERSAL_PAIRS = ((0, 5), (1, 4), (2, 3))
METHODS = (
    "crossfit_oracle_bq",
    "active_solvai_bq",
    "uniform_solvai_bq",
    "fixed_solvai_bq",
    "generic_bq",
    "uniform_direct",
    "fixed_direct",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_inputs(responses: pd.DataFrame) -> None:
    paths = {
        "responses": PHASE2 / "dense_responses_prospective.parquet",
        "lock": ACTIVE / "release/DENSE_SENTINEL_CALIBRATION_LOCK.json",
        "config": ACTIVE / "configs/dense_sentinel_v1.json",
        "priors": ACTIVE / "results/phase1/phase1_response_predictions.parquet",
        "manifest": ACTIVE / "simulations/dense_pimd2/manifest.csv",
    }
    for key, path in paths.items():
        observed = sha256(path)
        if observed != EXPECTED_HASHES[key]:
            raise AssertionError(f"Frozen input hash mismatch for {key}: {observed}")
    records = (
        "\n".join(
            f"{row.molecule_id}\t{row['lambda']:.17g}\t{row.energy_sha256}"
            for _, row in responses.sort_values(["molecule_id", "lambda"]).iterrows()
        )
        + "\n"
    )
    aggregate = hashlib.sha256(records.encode()).hexdigest()
    if aggregate != EXPECTED_HASHES["energy_records"]:
        raise AssertionError(f"Raw-energy record hash mismatch: {aggregate}")
    for row in responses.to_dict("records"):
        if sha256(Path(str(row["energy_file"]))) != str(row["energy_sha256"]):
            raise AssertionError(f"Raw-energy hash mismatch: {row['energy_file']}")


def response_blocks(responses: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for row in responses.to_dict("records"):
        frame = read_energy(Path(str(row["energy_file"])))
        values = frame.dHdL.dropna().to_numpy(float)
        if len(values) != 51:
            raise AssertionError(f"Expected 51 frames: {row['energy_file']}")
        chunks = np.array_split(values, 4)
        for block_index, chunk in enumerate(chunks):
            records.append(
                {
                    "molecule_name": row["molecule_name"],
                    "molecule_id": row["molecule_id"],
                    "functional_group_family": row["functional_group_family"],
                    "lambda": float(row["lambda"]),
                    "block": block_index,
                    "first_frame": int(sum(len(part) for part in chunks[:block_index])),
                    "last_frame_exclusive": int(
                        sum(len(part) for part in chunks[: block_index + 1])
                    ),
                    "frames": len(chunk),
                    "mean_dhdl_kcal_mol": float(chunk.mean()),
                    "sd_dhdl_kcal_mol": float(chunk.std(ddof=1)),
                }
            )
    return pd.DataFrame(records)


def side_curve(
    blocks: pd.DataFrame,
    molecule: str,
    selected_blocks: tuple[int, int],
    grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    molecule_rows = blocks.loc[
        blocks.molecule_name.eq(molecule) & blocks.block.isin(selected_blocks)
    ]
    means: list[float] = []
    sems: list[float] = []
    for value in grid:
        rows = molecule_rows.loc[np.isclose(molecule_rows["lambda"], value)].sort_values("block")
        if len(rows) != 2:
            raise AssertionError((molecule, selected_blocks, value, len(rows)))
        # Block sizes differ only for the final block, so reconstruct the pooled
        # response exactly from block means and frame counts.
        pooled = float(np.average(rows.mean_dhdl_kcal_mol, weights=rows.frames))
        sem = float(rows.mean_dhdl_kcal_mol.to_numpy(float).std(ddof=1) / np.sqrt(2.0))
        means.append(pooled)
        sems.append(sem)
    return np.asarray(means), np.asarray(sems)


def load_priors(names: list[str], grid: np.ndarray) -> dict[str, np.ndarray]:
    frame = pd.read_parquet(ACTIVE / "results/phase1/phase1_response_predictions.parquet")
    frame = frame.loc[
        frame.partition.eq("standardized_exclusion_primary")
        & frame.repeat.eq(-1)
        & np.isclose(frame.trajectory_fraction, 1.0)
        & frame.component.eq("total")
        & frame.molecule_name.isin(names)
    ]
    priors: dict[str, np.ndarray] = {}
    for name, rows in frame.groupby("molecule_name"):
        rows = rows.sort_values("lambda")
        priors[name] = interpolate_three_point_prior(
            rows.predicted_structure_only.to_numpy(float), grid
        )
    if set(priors) != set(names):
        raise AssertionError("Missing frozen structure prior")
    return priors


def crossfit_oracle_order(
    training_curve: np.ndarray,
    training_sem: np.ndarray,
    prior: np.ndarray,
    covariance: np.ndarray,
    weights: np.ndarray,
    noise_multiplier: float,
) -> list[int]:
    chosen = list(INITIAL)
    remaining = set(range(len(training_curve))) - set(chosen)
    training_target = float(weights @ training_curve)
    order: list[int] = []
    while remaining:
        scores: list[tuple[float, int]] = []
        for candidate in remaining:
            indices = np.asarray(chosen + [candidate], dtype=int)
            posterior = condition_curve(
                prior,
                covariance,
                indices,
                training_curve[indices],
                noise_multiplier * training_sem[indices] ** 2,
            )
            estimate, _ = posterior.integral(weights)
            scores.append((abs(estimate - training_target), candidate))
        _, selected = min(scores)
        chosen.append(selected)
        order.append(selected)
        remaining.remove(selected)
    return order


def evaluate(
    *,
    grid: np.ndarray,
    weights: np.ndarray,
    evaluation_curve: np.ndarray,
    evaluation_sem: np.ndarray,
    prior: np.ndarray,
    covariance: np.ndarray,
    noise_multiplier: float,
    selected: list[int],
    direct: bool,
) -> dict[str, float]:
    indices = np.asarray(selected, dtype=int)
    target = -float(weights @ evaluation_curve)
    if direct:
        reconstructed = observed_pchip(grid, indices, evaluation_curve[indices])
        prediction = -float(weights @ reconstructed)
        posterior_sd = math.nan
        covered = math.nan
    else:
        posterior = condition_curve(
            prior,
            covariance,
            indices,
            evaluation_curve[indices],
            noise_multiplier * evaluation_sem[indices] ** 2,
        )
        estimate, posterior_sd = posterior.integral(weights)
        prediction = -estimate
        covered = float(abs(prediction - target) <= Z90 * posterior_sd)
    return {
        "dense_integral_kcal_mol": target,
        "predicted_integral_kcal_mol": prediction,
        "signed_error_kcal_mol": prediction - target,
        "absolute_error_kcal_mol": abs(prediction - target),
        "posterior_sd_kcal_mol": posterior_sd,
        "covered_90": covered,
    }


def bootstrap(values: np.ndarray, level: float) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    samples = values[rng.integers(0, len(values), size=(BOOTSTRAP_RESAMPLES, len(values)))].mean(
        axis=1
    )
    tail = (1.0 - level) / 2.0
    return tuple(float(value) for value in np.quantile(samples, [tail, 1.0 - tail]))


def jaccard(left: set[int], right: set[int]) -> float:
    return len(left & right) / len(left | right) if left | right else 1.0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    config = json.loads((ACTIVE / "configs/dense_sentinel_v1.json").read_text())
    lock = json.loads((ACTIVE / "release/DENSE_SENTINEL_CALIBRATION_LOCK.json").read_text())
    responses = pd.read_parquet(PHASE2 / "dense_responses_prospective.parquet")
    verify_inputs(responses)
    grid = np.asarray(config["lambda_grid"], dtype=float)
    weights = trapezoid_weights(grid)
    names = list(config["prospective_molecules"])
    priors = load_priors(names, grid)
    blocks = response_blocks(responses)
    blocks.to_parquet(OUT / "trajectory_block_responses.parquet", index=False)

    generic_prior = np.asarray(lock["generic_prior_mean"], dtype=float)
    generic_settings = lock["selection"]["generic"]
    solvai_settings = lock["selection"]["solvai"]
    expected_sem = np.asarray(lock["expected_sem_by_lambda"], dtype=float)
    generic_covariance = rbf_covariance(
        grid, generic_settings["amplitude"], generic_settings["lengthscale"]
    )
    generic_order = variance_reduction_order(
        generic_covariance,
        weights,
        INITIAL,
        expected_sem**2 * generic_settings["noise_inflation"],
    )

    records: list[dict[str, object]] = []
    schedule_records: list[dict[str, object]] = []
    dense_records: list[dict[str, object]] = []
    metadata = responses.drop_duplicates("molecule_name").set_index("molecule_name")
    for name in names:
        prior = priors[name]
        solvai_covariance = rbf_covariance(
            grid,
            solvai_settings["amplitude"],
            solvai_settings["lengthscale"],
            curvature_scale(prior, grid),
        )
        active_order = variance_reduction_order(
            solvai_covariance,
            weights,
            INITIAL,
            expected_sem**2 * solvai_settings["noise_inflation"],
        )
        uniform_order = maximin_order(grid, INITIAL)
        population_order = fixed_order(grid)
        for split_id, (selection_blocks, evaluation_blocks) in SPLITS.items():
            selection_curve, selection_sem = side_curve(blocks, name, selection_blocks, grid)
            evaluation_curve, evaluation_sem = side_curve(blocks, name, evaluation_blocks, grid)
            selection_integral = -float(weights @ selection_curve)
            evaluation_integral = -float(weights @ evaluation_curve)
            dense_records.append(
                {
                    "molecule_name": name,
                    "molecule_id": metadata.loc[name, "molecule_id"],
                    "split": split_id,
                    "selection_blocks": ",".join(map(str, selection_blocks)),
                    "evaluation_blocks": ",".join(map(str, evaluation_blocks)),
                    "selection_dense_integral_kcal_mol": selection_integral,
                    "evaluation_dense_integral_kcal_mol": evaluation_integral,
                    "absolute_dense_difference_kcal_mol": abs(
                        evaluation_integral - selection_integral
                    ),
                }
            )
            oracle_order = crossfit_oracle_order(
                selection_curve,
                selection_sem,
                prior,
                solvai_covariance,
                weights,
                solvai_settings["noise_inflation"],
            )
            schedules = {
                "crossfit_oracle_bq": (
                    oracle_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    False,
                ),
                "active_solvai_bq": (
                    active_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    False,
                ),
                "uniform_solvai_bq": (
                    uniform_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    False,
                ),
                "fixed_solvai_bq": (
                    population_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    False,
                ),
                "generic_bq": (
                    generic_order,
                    generic_prior,
                    generic_covariance,
                    generic_settings["noise_inflation"],
                    False,
                ),
                "uniform_direct": (
                    uniform_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    True,
                ),
                "fixed_direct": (
                    population_order,
                    prior,
                    solvai_covariance,
                    solvai_settings["noise_inflation"],
                    True,
                ),
            }
            for method, (order, method_prior, kernel, multiplier, direct) in schedules.items():
                for budget in (5, 7):
                    selected = list(INITIAL) + order[: budget - len(INITIAL)]
                    result = evaluate(
                        grid=grid,
                        weights=weights,
                        evaluation_curve=evaluation_curve,
                        evaluation_sem=evaluation_sem,
                        prior=method_prior,
                        covariance=kernel,
                        noise_multiplier=multiplier,
                        selected=selected,
                        direct=direct,
                    )
                    records.append(
                        {
                            "molecule_name": name,
                            "molecule_id": metadata.loc[name, "molecule_id"],
                            "functional_group_family": metadata.loc[
                                name, "functional_group_family"
                            ],
                            "split": split_id,
                            "selection_blocks": ",".join(map(str, selection_blocks)),
                            "evaluation_blocks": ",".join(map(str, evaluation_blocks)),
                            "method": method,
                            "total_windows": budget,
                            "observed_lambdas": ",".join(f"{grid[index]:g}" for index in selected),
                            **result,
                        }
                    )
                    if method == "crossfit_oracle_bq":
                        schedule_records.append(
                            {
                                "molecule_name": name,
                                "molecule_id": metadata.loc[name, "molecule_id"],
                                "split": split_id,
                                "total_windows": budget,
                                "added_lambdas": ",".join(
                                    f"{grid[index]:g}" for index in selected[len(INITIAL) :]
                                ),
                                "added_indices": ",".join(
                                    str(index) for index in selected[len(INITIAL) :]
                                ),
                            }
                        )

    predictions = pd.DataFrame(records)
    schedules = pd.DataFrame(schedule_records)
    dense = pd.DataFrame(dense_records)
    predictions.to_parquet(OUT / "cross_block_predictions.parquet", index=False)
    schedules.to_csv(OUT / "crossfit_oracle_schedules.csv", index=False)
    dense.to_csv(OUT / "dense_integral_reproducibility_by_split.csv", index=False)

    molecule_metrics = predictions.groupby(
        ["method", "total_windows", "molecule_name", "molecule_id"], as_index=False
    ).agg(
        mean_absolute_error_kcal_mol=("absolute_error_kcal_mol", "mean"),
        mean_signed_error_kcal_mol=("signed_error_kcal_mol", "mean"),
        split_sd_absolute_error_kcal_mol=("absolute_error_kcal_mol", "std"),
        coverage_90=("covered_90", "mean"),
    )
    molecule_metrics.to_csv(OUT / "molecule_level_metrics.csv", index=False)

    aggregate_rows: list[dict[str, object]] = []
    for (method, budget), rows in molecule_metrics.groupby(["method", "total_windows"]):
        values = rows.mean_absolute_error_kcal_mol.to_numpy(float)
        low90, high90 = bootstrap(values, 0.90)
        low95, high95 = bootstrap(values, 0.95)
        aggregate_rows.append(
            {
                "method": method,
                "total_windows": budget,
                "n_molecules": len(rows),
                "mae_kcal_mol": float(values.mean()),
                "ci90_low": low90,
                "ci90_high": high90,
                "ci95_low": low95,
                "ci95_high": high95,
                "mean_coverage_90": float(rows.coverage_90.mean()),
            }
        )
    aggregate = pd.DataFrame(aggregate_rows)
    aggregate.to_csv(OUT / "aggregate_metrics.csv", index=False)

    comparison_rows: list[dict[str, object]] = []
    for budget in (5, 7):
        oracle = molecule_metrics.loc[
            molecule_metrics.method.eq("crossfit_oracle_bq")
            & molecule_metrics.total_windows.eq(budget)
        ].set_index("molecule_id")
        for comparator in METHODS[1:]:
            control = molecule_metrics.loc[
                molecule_metrics.method.eq(comparator) & molecule_metrics.total_windows.eq(budget)
            ].set_index("molecule_id")
            difference = (
                oracle.mean_absolute_error_kcal_mol - control.mean_absolute_error_kcal_mol
            ).to_numpy(float)
            low90, high90 = bootstrap(difference, 0.90)
            low95, high95 = bootstrap(difference, 0.95)
            comparison_rows.append(
                {
                    "total_windows": budget,
                    "candidate": "crossfit_oracle_bq",
                    "comparator": comparator,
                    "candidate_minus_comparator_mae_kcal_mol": float(difference.mean()),
                    "ci90_low": low90,
                    "ci90_high": high90,
                    "ci95_low": low95,
                    "ci95_high": high95,
                    "fraction_molecules_improved": float(np.mean(difference < 0)),
                }
            )
    comparisons = pd.DataFrame(comparison_rows)
    comparisons.to_csv(OUT / "paired_comparisons.csv", index=False)

    unique_dense = dense.loc[dense.split.isin([0, 1, 2])].copy()
    dense_molecule = unique_dense.groupby(
        ["molecule_name", "molecule_id"], as_index=False
    ).absolute_dense_difference_kcal_mol.mean()
    dense_values = dense_molecule.absolute_dense_difference_kcal_mol.to_numpy(float)
    dense_low90, dense_high90 = bootstrap(dense_values, 0.90)
    dense_low95, dense_high95 = bootstrap(dense_values, 0.95)
    squared_by_molecule = unique_dense.groupby("molecule_id").apply(
        lambda rows: float(np.mean(rows.absolute_dense_difference_kcal_mol.to_numpy(float) ** 2)),
        include_groups=False,
    )
    dense_reproducibility = {
        "n_molecules": len(dense_molecule),
        "complementary_partitions_per_molecule": 3,
        "mae_kcal_mol": float(dense_values.mean()),
        "median_absolute_difference_kcal_mol": float(
            unique_dense.absolute_dense_difference_kcal_mol.median()
        ),
        "rmse_kcal_mol": float(np.sqrt(squared_by_molecule.mean())),
        "ci90_low": dense_low90,
        "ci90_high": dense_high90,
        "ci95_low": dense_low95,
        "ci95_high": dense_high95,
    }

    jaccard_rows: list[dict[str, object]] = []
    for budget in (5, 7):
        budget_schedules = schedules.loc[schedules.total_windows.eq(budget)]
        for name in names:
            rows = budget_schedules.loc[budget_schedules.molecule_name.eq(name)].set_index("split")
            sets = {
                int(split): {int(value) for value in str(row.added_indices).split(",")}
                for split, row in rows.iterrows()
            }
            for left, right in REVERSAL_PAIRS:
                jaccard_rows.append(
                    {
                        "molecule_name": name,
                        "total_windows": budget,
                        "comparison": "reversal",
                        "split_left": left,
                        "split_right": right,
                        "jaccard_added_lambdas": jaccard(sets[left], sets[right]),
                    }
                )
            for left, right in itertools.combinations(range(6), 2):
                jaccard_rows.append(
                    {
                        "molecule_name": name,
                        "total_windows": budget,
                        "comparison": "all_pairs",
                        "split_left": left,
                        "split_right": right,
                        "jaccard_added_lambdas": jaccard(sets[left], sets[right]),
                    }
                )
    jaccards = pd.DataFrame(jaccard_rows)
    jaccards.to_csv(OUT / "oracle_schedule_stability.csv", index=False)

    original = pd.read_csv(PHASE2 / "dense_replay_metrics.csv")
    decision_rows: list[dict[str, object]] = []
    b_passes: list[bool] = []
    a_passes: list[bool] = []
    for budget in (5, 7):
        metric_lookup = aggregate.loc[aggregate.total_windows.eq(budget)].set_index("method")
        comparison = comparisons.loc[
            comparisons.total_windows.eq(budget) & comparisons.comparator.eq("uniform_direct")
        ].iloc[0]
        original_lookup = original.loc[original.total_windows.eq(budget)].set_index("method")
        original_headroom = float(
            original_lookup.loc["uniform_direct", "integral_mae_kcal_mol"]
            - original_lookup.loc["oracle_non_deployable", "integral_mae_kcal_mol"]
        )
        crossfit_headroom = float(
            metric_lookup.loc["uniform_direct", "mae_kcal_mol"]
            - metric_lookup.loc["crossfit_oracle_bq", "mae_kcal_mol"]
        )
        survival = crossfit_headroom / original_headroom
        b_pass = bool(
            comparison.candidate_minus_comparator_mae_kcal_mol <= -0.10
            and comparison.ci90_high < 0.0
            and comparison.fraction_molecules_improved >= 0.75
            and survival >= 0.50
        )
        a_pass = bool(
            comparison.candidate_minus_comparator_mae_kcal_mol >= -0.10
            and comparison.ci90_low > -0.10
            and survival <= 0.25
        )
        b_passes.append(b_pass)
        a_passes.append(a_pass)
        decision_rows.append(
            {
                "total_windows": budget,
                "crossfit_oracle_mae_kcal_mol": float(
                    metric_lookup.loc["crossfit_oracle_bq", "mae_kcal_mol"]
                ),
                "uniform_direct_mae_kcal_mol": float(
                    metric_lookup.loc["uniform_direct", "mae_kcal_mol"]
                ),
                "oracle_minus_uniform_kcal_mol": float(
                    comparison.candidate_minus_comparator_mae_kcal_mol
                ),
                "ci90_low": float(comparison.ci90_low),
                "ci90_high": float(comparison.ci90_high),
                "fraction_molecules_improved": float(comparison.fraction_molecules_improved),
                "original_same_curve_headroom_kcal_mol": original_headroom,
                "crossfit_headroom_kcal_mol": crossfit_headroom,
                "headroom_survival_ratio": survival,
                "a_budget_condition": a_pass,
                "b_budget_condition": b_pass,
            }
        )
    j7 = float(
        jaccards.loc[
            jaccards.total_windows.eq(7) & jaccards.comparison.eq("reversal"),
            "jaccard_added_lambdas",
        ].mean()
    )
    b_final = all(b_passes) and j7 >= 0.40
    a_final = all(a_passes)
    noisy_reference = bool(
        dense_reproducibility["mae_kcal_mol"] > 0.50
        or dense_reproducibility["median_absolute_difference_kcal_mol"] > 0.30
    )
    if b_final:
        conclusion = "B. STABLE BUT UNLEARNED HEADROOM"
    elif a_final:
        conclusion = "A. NO-GO CONFIRMED"
    else:
        conclusion = "C. INCONCLUSIVE DUE TO RESPONSE NOISE"

    metrics_payload = {
        "schema_version": 1,
        "protocol_commit": "d0167dcdfc4298a3b0a0ffbfa0be1b87bc8f2be3",
        "conclusion": conclusion,
        "dense_integral_reproducibility": dense_reproducibility,
        "decision_rows": decision_rows,
        "mean_reversal_jaccard_budget5": float(
            jaccards.loc[
                jaccards.total_windows.eq(5) & jaccards.comparison.eq("reversal"),
                "jaccard_added_lambdas",
            ].mean()
        ),
        "mean_reversal_jaccard_budget7": j7,
        "mean_all_pairs_jaccard_budget5": float(
            jaccards.loc[
                jaccards.total_windows.eq(5) & jaccards.comparison.eq("all_pairs"),
                "jaccard_added_lambdas",
            ].mean()
        ),
        "mean_all_pairs_jaccard_budget7": float(
            jaccards.loc[
                jaccards.total_windows.eq(7) & jaccards.comparison.eq("all_pairs"),
                "jaccard_added_lambdas",
            ].mean()
        ),
        "reliability_gate_exceeded": noisy_reference,
        "b_rule_passed": b_final,
        "a_rule_passed": a_final,
        "input_hashes": EXPECTED_HASHES,
    }
    metrics_path = OUT / "canonical_metrics.json"
    metrics_path.write_text(json.dumps(metrics_payload, indent=2) + "\n")

    # Diagnostic figures use only frozen outputs and the predeclared aggregations.
    colors = {
        "crossfit_oracle_bq": "#D55E00",
        "active_solvai_bq": "#0072B2",
        "uniform_direct": "#009E73",
        "generic_bq": "#777777",
        "fixed_direct": "#E69F00",
    }
    figure_methods = list(colors)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    for axis, budget in zip(axes, (5, 7), strict=True):
        plot = molecule_metrics.loc[
            molecule_metrics.total_windows.eq(budget) & molecule_metrics.method.isin(figure_methods)
        ]
        pivot = plot.pivot(
            index="molecule_name", columns="method", values="mean_absolute_error_kcal_mol"
        )
        x = np.arange(len(pivot))
        width = 0.15
        for offset, method in zip(
            np.linspace(-2 * width, 2 * width, len(figure_methods)),
            figure_methods,
            strict=True,
        ):
            axis.bar(x + offset, pivot[method], width, color=colors[method], label=method)
        axis.set_title(f"{budget} windows")
        axis.set_xticks(x, pivot.index, rotation=45, ha="right", fontsize=7)
        axis.set_ylabel("Held-out-block integral error (kcal mol$^{-1}$)")
        axis.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center", fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.80))
    for suffix in ("svg", "png"):
        fig.savefig(FIGURES / f"oracle_independent_noise_molecule_errors.{suffix}", dpi=300)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(5.7, 3.8), constrained_layout=True)
    decision = pd.DataFrame(decision_rows)
    x = np.arange(2)
    original_headroom = decision.original_same_curve_headroom_kcal_mol.to_numpy(float)
    crossfit_headroom = decision.crossfit_headroom_kcal_mol.to_numpy(float)
    axis.bar(x - 0.18, original_headroom, 0.36, color="#CC79A7", label="same-curve oracle")
    axis.bar(x + 0.18, crossfit_headroom, 0.36, color="#D55E00", label="cross-fitted oracle")
    axis.axhline(0, color="#333333", linewidth=0.8)
    axis.set_xticks(x, ["5 windows", "7 windows"])
    axis.set_ylabel("Headroom vs uniform direct (kcal mol$^{-1}$)")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    for suffix in ("svg", "png"):
        fig.savefig(FIGURES / f"oracle_headroom_survival.{suffix}", dpi=300)
    plt.close(fig)

    output_paths = sorted(OUT.glob("*")) + sorted(FIGURES.glob("oracle_*")) + [REPORT]
    artifact_manifest = {
        str(path.relative_to(ROOT)): {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in output_paths
        if path.is_file() and path.name != "artifact_manifest.json"
    }
    manifest_path = OUT / "artifact_manifest.json"
    manifest_path.write_text(json.dumps(artifact_manifest, indent=2) + "\n")
    ledger_path = ACTIVE / "runs/ledger.jsonl"
    run_id = "AS-V2-ORACLE-XBLOCK-001"
    existing_ids = {
        json.loads(line).get("run_id")
        for line in ledger_path.read_text().splitlines()
        if line.strip()
    }
    if run_id not in existing_ids:
        append_record(
            ledger_path,
            {
                "run_id": run_id,
                "timestamp_utc": datetime.now(UTC).isoformat(),
                "stage": "v2_oracle_independent_noise_diagnostic",
                "action": "six-way non-overlapping cross-block oracle audit",
                "status": "completed",
                "quality_control": conclusion,
                "device": "CPU",
                "gpu_hours": 0.0,
                "cpu_hours": 0.0,
                "simulated_time_ps": 0.0,
                "bead_windows": 0,
                "force_evaluations": 0,
                "protocol_commit": "d0167dcdfc4298a3b0a0ffbfa0be1b87bc8f2be3",
                "output": str(metrics_path),
                "output_sha256": sha256(metrics_path),
            },
        )

    print(json.dumps(metrics_payload, indent=2))


if __name__ == "__main__":
    main()
