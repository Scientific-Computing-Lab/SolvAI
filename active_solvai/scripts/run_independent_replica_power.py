"""Prospective power gate for the frozen independent-replica experiment.

This script implements only the calculation registered in
INDEPENDENT_REPLICA_RESOLUTION_FREEZE.md.  It cannot launch simulations.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai/results/v2_independent_replicas/power"
REPORT = ROOT / "active_solvai/reports/INDEPENDENT_REPLICA_POWER_ANALYSIS.md"
SEED = 20260904
N_TRIALS = 100_000
TARGET_EFFECT = -0.30
THROUGHPUT_PS_GPU_H = 369.902

INPUTS = {
    "cross_block_predictions": (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/cross_block_predictions.parquet",
        "e6492f5754aab9d4724e10c9b5e2d08863e669beb9b252a4f68a3e0b88a715fa",
    ),
    "molecule_level_metrics": (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/molecule_level_metrics.csv",
        "7be6b7c3580b8386e6f31395463a9a6d6723d41c1590f2d714a066dd4fff3679",
    ),
    "dense_reproducibility": (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/dense_integral_reproducibility_by_split.csv",
        "feec7765a4705ca9390e87fd83eeb29847539cbc18f2aaa16f08bc7c71c3df0e",
    ),
    "schedule_stability": (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/oracle_schedule_stability.csv",
        "09a30eb8afacd2220544e4f6f07163e8655fa2b6975e3a2501669d80eda2cdeb",
    ),
    "canonical_metrics": (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/canonical_metrics.json",
        "c0f01f65f3cfae4230d01fb536b7f375a1a0db29742a5023870ed53022cec2f4",
    ),
    "responses": (
        ROOT / "active_solvai/results/phase2/dense_responses_prospective.parquet",
        "4ce191a2ebfabb10dfcf5e11f98fef91b0bd5b2993f066ab47cfa9cc261c1097",
    ),
    "replay": (
        ROOT / "active_solvai/results/phase2/dense_replay_predictions.parquet",
        "a0dbc9b61c7f54848be8e56dcc5c22d12d511839c90cdbfe44900574112bf1fd",
    ),
    "config": (
        ROOT / "active_solvai/configs/dense_sentinel_v1.json",
        "fc6afd2a4f8a255a641228df624769e6ed437dc4cc69dc244480a33c2936594d",
    ),
    "manifest": (
        ROOT / "active_solvai/simulations/dense_pimd2/manifest.csv",
        "9c46ad5b78af027f3c82d84ff300c1ea839fbd13650eba704bbfef95f5871358",
    ),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_inputs() -> dict[str, str]:
    observed = {}
    for key, (path, expected) in INPUTS.items():
        digest = sha256(path)
        if digest != expected:
            raise RuntimeError(f"hash mismatch for {key}: {digest} != {expected}")
        observed[key] = digest
    return observed


def effect_components() -> tuple[dict[int, np.ndarray], dict[int, float], pd.DataFrame]:
    frame = pd.read_parquet(INPUTS["cross_block_predictions"][0])
    keep = frame[frame["method"].isin(["crossfit_oracle_bq", "uniform_direct"])].copy()
    pivot = keep.pivot_table(
        index=["molecule_name", "molecule_id", "split", "total_windows"],
        columns="method",
        values="absolute_error_kcal_mol",
        aggfunc="first",
    ).reset_index()
    pivot["effect"] = pivot["crossfit_oracle_bq"] - pivot["uniform_direct"]
    reversal_partition = {0: 0, 5: 0, 1: 1, 4: 1, 2: 2, 3: 2}
    pivot["partition"] = pivot["split"].map(reversal_partition)
    paired = (
        pivot.groupby(
            ["molecule_name", "molecule_id", "total_windows", "partition"], as_index=False
        )["effect"]
        .mean()
        .sort_values(["total_windows", "molecule_name", "partition"])
    )

    molecule_patterns: dict[int, np.ndarray] = {}
    residual_variances: dict[int, float] = {}
    component_rows = []
    for budget in (5, 7):
        b = paired[paired.total_windows == budget]
        matrix = b.pivot(index="molecule_name", columns="partition", values="effect").sort_index()
        means = matrix.mean(axis=1).to_numpy(float)
        centered = means - means.mean()
        residuals = matrix.to_numpy(float) - means[:, None]
        residual_variance = float(np.sum(residuals**2) / (residuals.size - len(means)))
        molecule_patterns[budget] = centered
        residual_variances[budget] = residual_variance
        for molecule, mean, centered_value, row in zip(
            matrix.index, means, centered, matrix.to_numpy(float), strict=True
        ):
            component_rows.append(
                {
                    "budget": budget,
                    "molecule_name": molecule,
                    "observed_mean_effect": mean,
                    "centered_effect": centered_value,
                    "target_shifted_effect": centered_value + TARGET_EFFECT,
                    "partition_0_effect": row[0],
                    "partition_1_effect": row[1],
                    "partition_2_effect": row[2],
                    "pooled_residual_variance_2p5ps": residual_variance,
                }
            )
    return molecule_patterns, residual_variances, pd.DataFrame(component_rows)


def dense_scales() -> tuple[np.ndarray, pd.DataFrame]:
    frame = pd.read_csv(INPUTS["dense_reproducibility"][0])
    unique = frame[frame["split"].isin([0, 1, 2])].copy()
    matrix = unique.pivot(
        index="molecule_name", columns="split", values="absolute_dense_difference_kcal_mol"
    ).sort_index()
    scales = matrix.mean(axis=1).to_numpy(float)
    detail = matrix.reset_index()
    detail["mean_abs_dense_difference_2p5ps"] = scales
    return scales, detail


def draw_patterns(
    rng: np.random.Generator, base: np.ndarray, n_molecules: int
) -> np.ndarray:
    if n_molecules == 8:
        return np.broadcast_to(base, (N_TRIALS, 8)).copy()
    indices = rng.integers(0, len(base), size=(N_TRIALS, n_molecules))
    return base[indices]


def power_for_design(
    rng: np.random.Generator,
    centered: np.ndarray,
    residual_variance: float,
    n_molecules: int,
    production_ps: float,
    variance_multiplier: float = 1.0,
) -> float:
    true_effect = draw_patterns(rng, centered, n_molecules) + TARGET_EFFECT
    direction_sd = math.sqrt(residual_variance * (2.5 / production_ps) * variance_multiplier)
    # The registered endpoint averages the two independent transfer directions.
    observed = true_effect + rng.normal(
        0.0, direction_sd / math.sqrt(2.0), size=true_effect.shape
    )
    means = observed.mean(axis=1)
    sds = observed.std(axis=1, ddof=1)
    upper90 = means + stats.t.ppf(0.90, n_molecules - 1) * sds / math.sqrt(n_molecules)
    fraction_improved = (observed < 0).mean(axis=1)
    detected = (means <= -0.20) & (upper90 < 0.0) & (fraction_improved >= 0.75)
    return float(detected.mean())


def mean_only_power_for_design(
    rng: np.random.Generator,
    centered: np.ndarray,
    residual_variance: float,
    n_molecules: int,
    production_ps: float,
) -> float:
    """Planning-only power for a mean effect, without the 75% consistency gate."""
    true_effect = draw_patterns(rng, centered, n_molecules) + TARGET_EFFECT
    direction_sd = math.sqrt(residual_variance * (2.5 / production_ps))
    observed = true_effect + rng.normal(
        0.0, direction_sd / math.sqrt(2.0), size=true_effect.shape
    )
    means = observed.mean(axis=1)
    sds = observed.std(axis=1, ddof=1)
    upper90 = means + stats.t.ppf(0.90, n_molecules - 1) * sds / math.sqrt(n_molecules)
    return float(((means <= -0.20) & (upper90 < 0.0)).mean())


def reliability_for_design(
    rng: np.random.Generator,
    base_scales: np.ndarray,
    n_molecules: int,
    production_ps: float,
) -> float:
    if n_molecules == 8:
        scales = np.broadcast_to(base_scales, (N_TRIALS, 8)).copy()
    else:
        idx = rng.integers(0, len(base_scales), size=(N_TRIALS, n_molecules))
        scales = base_scales[idx]
    # Normalize |N(0,1)| to mean one, so each molecule retains its observed
    # complementary-half absolute-difference scale before inverse-time scaling.
    half_normal = np.abs(rng.normal(size=scales.shape)) / math.sqrt(2.0 / math.pi)
    differences = scales * math.sqrt(2.5 / production_ps) * half_normal
    passed = (differences.mean(axis=1) <= 0.50) & (
        np.median(differences, axis=1) <= 0.30
    )
    return float(passed.mean())


def evaluate_grid(
    patterns: dict[int, np.ndarray], residual_variances: dict[int, float], scales: np.ndarray
) -> pd.DataFrame:
    rows = []
    for n_molecules in (8, 12, 16, 24, 32, 48, 64, 96):
        for production_ps in (50, 75, 100, 150, 200):
            # Seed every design/budget independently so the output does not depend
            # on loop execution details.
            p = {}
            pc = {}
            for budget in (5, 7):
                rng = np.random.default_rng(
                    SEED + n_molecules * 10000 + production_ps * 10 + budget
                )
                p[budget] = power_for_design(
                    rng,
                    patterns[budget],
                    residual_variances[budget],
                    n_molecules,
                    production_ps,
                    variance_multiplier=1.0,
                )
                rng_c = np.random.default_rng(
                    SEED + 5_000_000 + n_molecules * 10000 + production_ps * 10 + budget
                )
                pc[budget] = power_for_design(
                    rng_c,
                    patterns[budget],
                    residual_variances[budget],
                    n_molecules,
                    production_ps,
                    variance_multiplier=2.0,
                )
            rng_rel = np.random.default_rng(
                SEED + 9_000_000 + n_molecules * 10000 + production_ps * 10
            )
            reliability = reliability_for_design(rng_rel, scales, n_molecules, production_ps)
            production_total_ps = n_molecules * 15 * 2 * production_ps
            gpu_hours = production_total_ps / THROUGHPUT_PS_GPU_H
            rows.append(
                {
                    "n_molecules": n_molecules,
                    "production_ps_per_replica": production_ps,
                    "replicas": 2,
                    "lambda_windows": 15,
                    "target_mean_effect_kcal_mol": TARGET_EFFECT,
                    "power_budget_5": p[5],
                    "power_budget_7": p[7],
                    "conservative_power_budget_5": pc[5],
                    "conservative_power_budget_7": pc[7],
                    "dense_reliability_probability": reliability,
                    "adequately_powered": bool(
                        p[5] >= 0.80 and p[7] >= 0.80 and reliability >= 0.80
                    ),
                    "production_total_ps": production_total_ps,
                    "projected_production_gpu_hours": gpu_hours,
                    "operational_reservation_gpu_hours": gpu_hours * 1.11,
                }
            )
    return pd.DataFrame(rows)


def evaluate_mean_only_planning_boundary(
    patterns: dict[int, np.ndarray], residual_variances: dict[int, float], scales: np.ndarray
) -> pd.DataFrame:
    """Quantify the least-cost mean-only resolution boundary.

    This does not amend the registered launch criterion.  The candidates cover
    the first passing 75-ps design and the maximum molecule counts at every
    other duration that cost less than that design.
    """
    candidates = [(157, 50), (104, 75), (105, 75), (78, 100), (52, 150), (39, 200)]
    rows = []
    for n_molecules, production_ps in candidates:
        power = {}
        for budget in (5, 7):
            rng = np.random.default_rng(
                SEED + 20_000_000 + n_molecules * 10000 + production_ps * 10 + budget
            )
            power[budget] = mean_only_power_for_design(
                rng,
                patterns[budget],
                residual_variances[budget],
                n_molecules,
                production_ps,
            )
        rng_rel = np.random.default_rng(
            SEED + 30_000_000 + n_molecules * 10000 + production_ps * 10
        )
        reliability = reliability_for_design(rng_rel, scales, n_molecules, production_ps)
        gpu_hours = n_molecules * 15 * 2 * production_ps / THROUGHPUT_PS_GPU_H
        rows.append(
            {
                "n_molecules": n_molecules,
                "production_ps_per_replica": production_ps,
                "mean_only_power_budget_5": power[5],
                "mean_only_power_budget_7": power[7],
                "dense_reliability_probability": reliability,
                "mean_only_adequately_powered": bool(
                    power[5] >= 0.80 and power[7] >= 0.80 and reliability >= 0.80
                ),
                "projected_production_gpu_hours": gpu_hours,
                "operational_reservation_gpu_hours": gpu_hours * 1.11,
                "supports_registered_consistency_claim": False,
            }
        )
    return pd.DataFrame(rows)


def make_power_figure(grid: pd.DataFrame, components: pd.DataFrame) -> None:
    figure_dir = ROOT / "active_solvai/figures/v2_independent_replicas"
    figure_dir.mkdir(parents=True, exist_ok=True)
    eight = grid[grid.n_molecules == 8].sort_values("production_ps_per_replica")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)

    ax = axes[0]
    ax.plot(
        eight.production_ps_per_replica,
        eight.power_budget_5,
        marker="o",
        color="#0072B2",
        label="Placement power, 5 windows",
    )
    ax.plot(
        eight.production_ps_per_replica,
        eight.power_budget_7,
        marker="s",
        color="#009E73",
        label="Placement power, 7 windows",
    )
    ax.plot(
        eight.production_ps_per_replica,
        eight.dense_reliability_probability,
        marker="^",
        color="#D55E00",
        label="Dense-reference reliability",
    )
    ax.axhline(0.8, color="#555555", lw=1.0, ls="--", label="Registered gate")
    ax.scatter([50], [eight.iloc[0].dense_reliability_probability], s=75, facecolors="none", edgecolors="#222222", zorder=5)
    ax.set(xlabel="Production per independent replica (ps)", ylabel="Probability", ylim=(-0.03, 1.03), title="a  Eight-sentinel precision gate")
    ax.legend(frameon=False, fontsize=8, loc="center right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    b5 = components[components.budget == 5].sort_values("target_shifted_effect")
    order = b5.molecule_name.tolist()
    for j, budget in enumerate((5, 7)):
        values = (
            components[components.budget == budget]
            .set_index("molecule_name")
            .loc[order, "target_shifted_effect"]
            .to_numpy()
        )
        y = np.arange(len(order)) + (j - 0.5) * 0.22
        ax.scatter(
            values,
            y,
            s=34,
            marker="o" if budget == 5 else "s",
            color="#0072B2" if budget == 5 else "#009E73",
            label=f"{budget} windows",
        )
    ax.axvline(0, color="#555555", lw=1.0)
    ax.axvline(-0.3, color="#CC79A7", lw=1.0, ls="--", label="Target mean")
    ax.set_yticks(np.arange(len(order)), order)
    ax.set(xlabel="Oracle minus uniform absolute error (kcal mol$^{-1}$)", title="b  Empirical heterogeneity retained")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Independent-replica experiment is blocked before simulation", fontsize=12, fontweight="bold")
    for suffix in ("svg", "png"):
        fig.savefig(figure_dir / f"power_gate.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    hashes = assert_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    patterns, residual_variances, components = effect_components()
    scales, dense_detail = dense_scales()
    grid = evaluate_grid(patterns, residual_variances, scales)
    mean_only_boundary = evaluate_mean_only_planning_boundary(
        patterns, residual_variances, scales
    )
    proposed = grid[(grid.n_molecules == 8) & (grid.production_ps_per_replica == 50)].iloc[0]
    passing = grid[grid.adequately_powered].sort_values(
        ["projected_production_gpu_hours", "n_molecules", "production_ps_per_replica"]
    )
    minimum = None if passing.empty else passing.iloc[0]
    launch = bool(proposed.adequately_powered)

    components.to_csv(OUT / "effect_variance_components.csv", index=False)
    dense_detail.to_csv(OUT / "dense_reliability_inputs.csv", index=False)
    grid.to_csv(OUT / "power_design_grid.csv", index=False)
    mean_only_boundary.to_csv(OUT / "mean_only_planning_boundary.csv", index=False)
    make_power_figure(grid, components)

    payload = {
        "schema_version": 1,
        "protocol_commit": "52827f851ed0d1a540b0492c63c3df269cf774e1",
        "monte_carlo_seed": SEED,
        "monte_carlo_trials": N_TRIALS,
        "target_mean_effect_kcal_mol": TARGET_EFFECT,
        "input_hashes": hashes,
        "pooled_residual_variance_2p5ps": {
            str(k): v for k, v in residual_variances.items()
        },
        "proposed_design": proposed.to_dict(),
        "launch_authorized_by_gate": launch,
        "smallest_adequately_powered_alternative": (
            None if minimum is None else minimum.to_dict()
        ),
        "planning_only_mean_effect_alternative": mean_only_boundary[
            mean_only_boundary.mean_only_adequately_powered
        ]
        .sort_values("projected_production_gpu_hours")
        .iloc[0]
        .to_dict(),
        "simulation_started": False,
    }
    (OUT / "power_analysis.json").write_text(json.dumps(payload, indent=2, default=float) + "\n")

    if minimum is None:
        alt_text = "No design on the frozen grid passed all power and reliability criteria."
    else:
        alt_text = (
            f"The least-cost passing grid design is **{int(minimum.n_molecules)} molecules, "
            f"two replicas, 15 windows and {int(minimum.production_ps_per_replica)} ps per replica**: "
            f"{minimum.projected_production_gpu_hours:.1f} projected production GPU-h and "
            f"{minimum.operational_reservation_gpu_hours:.1f} reserved GPU-h."
        )
    decision = (
        "PASS — the registered eight-molecule simulation may launch."
        if launch
        else "FAIL — the registered eight-molecule simulation must not launch."
    )
    mean_alt = payload["planning_only_mean_effect_alternative"]
    report = f"""# Power and precision gate for independent PIMD2 replicas

**Prospective protocol:** commit `52827f851ed0d1a540b0492c63c3df269cf774e1`  
**Decision:** {decision}

This calculation precedes every new simulation. It asks whether two 50-ps replicas on the eight frozen sentinels can resolve a true cross-replica oracle advantage of 0.30 kcal mol-1 at both five and seven windows while also producing a reliable dense reference.

## Registered eight-molecule design

| Quantity | Result |
|---|---:|
| Power, 5 windows | {proposed.power_budget_5:.3f} |
| Power, 7 windows | {proposed.power_budget_7:.3f} |
| Conservative power, 5 windows | {proposed.conservative_power_budget_5:.3f} |
| Conservative power, 7 windows | {proposed.conservative_power_budget_7:.3f} |
| Probability dense-reference gate passes | {proposed.dense_reliability_probability:.3f} |
| Required probability for each primary gate | 0.800 |
| Projected production GPU-hours | {proposed.projected_production_gpu_hours:.2f} |
| Operational reservation GPU-hours | {proposed.operational_reservation_gpu_hours:.2f} |

The design is launchable only when both budget-specific powers and dense reliability are at least 0.80. The decision is therefore mechanical rather than result-dependent.

## Smallest adequately powered alternative

{alt_text}

This alternative is **not authorized or launched**. If it requires more than eight molecules, those molecules must be chosen under a separate chemistry-first prospective freeze; this power calculation does not select them.

No amount of additional trajectory length or molecule replication on the registered empirical effect distribution makes the *full* decision rule adequately powered: after centering to a mean benefit of 0.30 kcal mol-1, only 4/8 molecules at five windows and 5/8 at seven windows retain favorable effects, below the required 75% consistency. Increased sample size therefore estimates that inconsistency more precisely rather than curing it.

For planning context only, dropping the consistency requirement and asking solely whether the **mean** effect is below zero gives a first passing boundary at **{int(mean_alt['n_molecules'])} molecules, two replicas, 15 windows and {int(mean_alt['production_ps_per_replica'])} ps per replica**. Its simulated powers are {mean_alt['mean_only_power_budget_5']:.3f} and {mean_alt['mean_only_power_budget_7']:.3f}, with dense-reference reliability probability {mean_alt['dense_reliability_probability']:.3f}. It requires {mean_alt['projected_production_gpu_hours']:.1f} projected production GPU-h and {mean_alt['operational_reservation_gpu_hours']:.1f} reserved GPU-h. This would resolve an average effect only; it cannot establish stable, broadly shared molecule-specific placement and is not recommended or authorized.

## Interpretation

The analysis preserves the observed between-molecule effect pattern, shifts only its mean to the scientifically meaningful alternative of -0.30 kcal mol-1, and projects finite-trajectory residual variance from 2.5 ps to the candidate duration by inverse-time scaling. The conservative sensitivity doubles residual variance. Dense-integral reliability is independently projected from the observed complementary-half discrepancies.

The calculation cannot prove that inverse-time scaling will hold or that a true oracle benefit exists. It only determines whether the proposed experiment has a reasonable chance to resolve an effect worth pursuing under the most favorable registered noise law. Failure therefore blocks simulation; passing merely authorizes the frozen experiment.

## Reproduction

```bash
uv run --project active_solvai python active_solvai/scripts/run_independent_replica_power.py
```

Machine-readable results are in `active_solvai/results/v2_independent_replicas/power/`.
"""
    REPORT.write_text(report)

    manifest = {}
    for path in sorted(OUT.glob("*")):
        if path.name == "artifact_manifest.json":
            continue
        manifest[str(path.relative_to(ROOT))] = sha256(path)
    manifest[str(REPORT.relative_to(ROOT))] = sha256(REPORT)
    for path in sorted((ROOT / "active_solvai/figures/v2_independent_replicas").glob("*")):
        manifest[str(path.relative_to(ROOT))] = sha256(path)
    (OUT / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(payload, indent=2, default=float))


if __name__ == "__main__":
    main()
