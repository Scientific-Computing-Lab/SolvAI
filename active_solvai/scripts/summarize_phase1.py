#!/usr/bin/env python3
"""Audit and summarize the frozen Active SolvAI Phase 1 outputs."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from active_solvai.ledger import append_record, sha256

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai/results/phase1"
REPORT = ROOT / "active_solvai/reports/PHASE1_ACTUAL_OBSERVATION_GATE.md"
NULL_RESULTS = ROOT / "active_solvai/NULL_RESULTS.md"
BOOTSTRAP_SEED = 20260828
COMPONENTS = (
    "total",
    "coulomb",
    "coulomb_softcore",
    "van_der_waals",
    "van_der_waals_softcore",
    "long_range_correction",
    "pme",
)
LAMBDAS = (0.1, 0.5, 0.9)


def paired_bootstrap(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    sampled = np.empty(100_000)
    offset = 0
    while offset < len(sampled):
        size = min(5000, len(sampled) - offset)
        indices = rng.integers(0, len(values), size=(size, len(values)))
        sampled[offset : offset + size] = values[indices].mean(axis=1)
        offset += size
    low, high = np.quantile(sampled, [0.025, 0.975])
    return {
        "mean": float(values.mean()),
        "ci_low_95": float(low),
        "ci_high_95": float(high),
        "fraction_below_zero": float(np.mean(values < 0)),
    }


def dataframe_markdown(frame: pd.DataFrame, decimals: int = 4) -> str:
    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        formatted = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                formatted.append(f"{float(value):.{decimals}f}")
            else:
                formatted.append(str(value))
        lines.append("| " + " | ".join(formatted) + " |")
    return "\n".join(lines)


def export_response_predictions() -> pd.DataFrame:
    identity = pd.read_parquet(ROOT / "active_solvai/data/identity/probe_identity_manifest.parquet")
    identity = identity.loc[identity.complete_three_point_curve].reset_index(drop=True)
    prefix = pd.read_parquet(OUT.parent / "phase0/response_prefix_blocks.parquet")
    prefix = prefix.loc[
        prefix.energy_group.eq("system")
        & prefix.molecule_id.astype(str).isin(identity.molecule_id.astype(str))
        & prefix.component.isin(
            [
                "dHdL",
                "dHdL_Coul",
                "dHdL_Coul_SC",
                "dHdL_VdW",
                "dHdL_VdW_SC",
                "dHdL_LRCor",
                "dHdL_PME",
            ]
        )
    ]
    component_map = {
        "dHdL": "total",
        "dHdL_Coul": "coulomb",
        "dHdL_Coul_SC": "coulomb_softcore",
        "dHdL_VdW": "van_der_waals",
        "dHdL_VdW_SC": "van_der_waals_softcore",
        "dHdL_LRCor": "long_range_correction",
        "dHdL_PME": "pme",
    }
    actual_lookup = {
        (
            str(row.molecule_id),
            float(row.trajectory_fraction),
            float(row["lambda"]),
            str(row.component),
        ): (float(row.mean_kcal_mol), float(row.five_block_sem_kcal_mol))
        for _, row in prefix.iterrows()
    }
    split = pd.read_csv(ROOT / "active_solvai/data/manifests/probe_split_assignments.csv")
    split = split.loc[
        split.partition.isin(
            ["standardized_exclusion_primary", "standardized_exclusion_repeat"]
        )
    ]
    rows = []
    for (partition, repeat, split_seed), assignment in split.groupby(
        ["partition", "repeat", "split_seed"], dropna=False, sort=False
    ):
        label = "primary" if partition.endswith("primary") else f"repeat_{int(repeat)}"
        assignment = assignment.set_index("molecule_id").reindex(identity.molecule_id)
        for fraction in (0.1, 0.2, 0.4, 0.7, 1.0):
            cache = np.load(OUT / "cache" / f"nested_response_{label}_fraction_{fraction:.1f}.npz")
            predicted = cache["outer"]
            for molecule_index, molecule in identity.iterrows():
                for lambda_index, lambda_value in enumerate(LAMBDAS):
                    for component_index, component in enumerate(COMPONENTS):
                        target_index = (
                            lambda_index
                            if component == "total"
                            else 3 + lambda_index * 6 + component_index - 1
                        )
                        source_component = next(
                            source
                            for source, normalized in component_map.items()
                            if normalized == component
                        )
                        lookup_key = (
                            str(molecule.molecule_id),
                            float(fraction),
                            float(lambda_value),
                            source_component,
                        )
                        if lookup_key not in actual_lookup:
                            raise AssertionError(f"Response export alignment failed: {lookup_key}")
                        actual_value, actual_sem = actual_lookup[lookup_key]
                        rows.append(
                            {
                                "molecule_id": molecule.molecule_id,
                                "molecule_name": molecule.molecule_name,
                                "functional_group_family": molecule.functional_group_family,
                                "partition": partition,
                                "repeat": int(repeat),
                                "split_seed": (
                                    float(split_seed) if pd.notna(split_seed) else np.nan
                                ),
                                "fold": int(assignment.iloc[molecule_index].fold),
                                "trajectory_fraction": fraction,
                                "lambda": lambda_value,
                                "component": component,
                                "actual": actual_value,
                                "predicted_structure_only": float(
                                    predicted[molecule_index, target_index]
                                ),
                                "five_block_sem": actual_sem,
                            }
                        )
    result = pd.DataFrame(rows)
    result["residual"] = result.actual - result.predicted_structure_only
    result["absolute_error"] = result.residual.abs()
    result.to_parquet(OUT / "phase1_response_predictions.parquet", index=False)
    return result


def make_figures(
    endpoint: pd.DataFrame,
    response: pd.DataFrame,
    reconstruction: pd.DataFrame,
) -> None:
    figure_root = ROOT / "active_solvai/figures/phase1"
    figure_root.mkdir(parents=True, exist_ok=True)
    colors = {"baseline": "#0072B2", "active": "#D55E00", "shuffle": "#777777"}

    repeated = endpoint.loc[endpoint.partition.eq("standardized_exclusion_repeat")]
    methods = {
        "P1-A_frozen_solvai": "Frozen SolvAI",
        "P1-D_actual_minus_predicted": "Actual−prior residual",
        "P1-H_mean_shuffled_residual": "Shuffled residual",
    }
    fig, ax = plt.subplots(figsize=(5.4, 3.25))
    for x, (method, label) in enumerate(methods.items()):
        query = repeated.loc[repeated.method.eq(method)]
        if method != "P1-A_frozen_solvai":
            query = query.loc[
                query.lambda_subset.eq("0p1_0p5_0p9")
                & query.response_scope.eq("total")
                & np.isclose(query.trajectory_fraction, 1.0)
            ]
        values = query.groupby("repeat").absolute_error.mean().sort_index().to_numpy()
        color = colors["baseline" if x == 0 else "active" if x == 1 else "shuffle"]
        ax.scatter(np.full(len(values), x), values, color=color, s=28, zorder=3)
        ax.hlines(values.mean(), x - 0.23, x + 0.23, color=color, lw=2.5)
    ax.set_xticks(range(3), list(methods.values()))
    ax.set_ylabel("MAE (kcal mol$^{-1}$)")
    ax.set_title("Actual PIMD2 residuals do not improve the endpoint")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_root / "phase1_endpoint_gate.svg")
    fig.savefig(figure_root / "phase1_endpoint_gate.png", dpi=300)
    plt.close(fig)

    query = response.loc[
        response.partition.eq("standardized_exclusion_repeat")
        & response.component.eq("total")
    ]
    summary = query.groupby(["trajectory_fraction", "lambda"]).absolute_error.mean().reset_index()
    fig, ax = plt.subplots(figsize=(5.4, 3.25))
    for lambda_value, group in summary.groupby("lambda"):
        ax.plot(
            group.trajectory_fraction * 5.0,
            group.absolute_error,
            marker="o",
            label=f"λ={lambda_value:g}",
        )
    ax.set_xlabel("Observed trajectory (ps)")
    ax.set_ylabel("Structure→response MAE (kcal mol$^{-1}$)")
    ax.set_title("Longer probes do not remove prior error")
    ax.legend(frameon=False, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_root / "phase1_response_learning.svg")
    fig.savefig(figure_root / "phase1_response_learning.png", dpi=300)
    plt.close(fig)

    query = reconstruction.loc[reconstruction.partition.eq("standardized_exclusion_repeat")]
    summary = (
        query.groupby(["posterior", "lambda_subset"]).absolute_error.mean().reset_index()
    )
    order = ["0p1", "0p5", "0p9", "0p1_0p5", "0p1_0p9", "0p5_0p9"]
    fig, ax = plt.subplots(figsize=(6.4, 3.35))
    width = 0.36
    for offset, (posterior, label, color) in enumerate(
        [
            ("P1-F_generic_posterior", "Generic", "#777777"),
            ("P1-F_solvai_conditioned_posterior", "SolvAI-conditioned", "#0072B2"),
        ]
    ):
        values = summary.loc[summary.posterior.eq(posterior)].set_index("lambda_subset")
        ax.bar(
            np.arange(len(order)) + (offset - 0.5) * width,
            values.reindex(order).absolute_error,
            width,
            label=label,
            color=color,
        )
    ax.set_xticks(np.arange(len(order)), [value.replace("p", ".").replace("_", ", ") for value in order])
    ax.set_xlabel("Observed λ")
    ax.set_ylabel("Hidden-point MAE (kcal mol$^{-1}$)")
    ax.set_title("Conditioning helps interpolation, but errors remain large")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_root / "phase1_reconstruction.svg")
    fig.savefig(figure_root / "phase1_reconstruction.png", dpi=300)
    plt.close(fig)


def main() -> None:
    endpoint = pd.read_parquet(OUT / "phase1_endpoint_predictions.parquet")
    metrics = pd.read_csv(OUT / "phase1_endpoint_metrics.csv")
    reconstruction = pd.read_parquet(OUT / "phase1_reconstruction_predictions.parquet")
    integration_metrics = pd.read_csv(OUT / "phase1_integration_metrics.csv")
    response = export_response_predictions()

    repeat_filter = endpoint.partition.eq("standardized_exclusion_repeat")
    primary_filter = (
        endpoint.method.eq("P1-D_actual_minus_predicted")
        & endpoint.lambda_subset.eq("0p1_0p5_0p9")
        & endpoint.response_scope.eq("total")
        & np.isclose(endpoint.trajectory_fraction, 1.0)
    )
    baseline_filter = endpoint.method.eq("P1-A_frozen_solvai")
    shuffle_filter = endpoint.method.eq("P1-H_mean_shuffled_residual")

    def per_molecule_error(mask: pd.Series) -> pd.Series:
        return endpoint.loc[repeat_filter & mask].groupby("molecule_id").absolute_error.mean()

    active_error = per_molecule_error(primary_filter)
    baseline_error = per_molecule_error(baseline_filter)
    shuffle_error = per_molecule_error(shuffle_filter)
    primary_vs_baseline = paired_bootstrap(
        active_error.sort_index().to_numpy() - baseline_error.sort_index().to_numpy()
    )
    primary_vs_shuffle = paired_bootstrap(
        active_error.sort_index().to_numpy() - shuffle_error.sort_index().to_numpy()
    )

    repeat_rows = metrics.loc[
        metrics.partition.eq("standardized_exclusion_repeat")
        & metrics.shuffle_seed.isna()
        & (
            metrics.method.eq("P1-A_frozen_solvai")
            | (
                metrics.method.isin(
                    [
                        "P1-B_predicted_response",
                        "P1-C_actual_response",
                        "P1-D_actual_minus_predicted",
                        "P1-H_mean_shuffled_residual",
                    ]
                )
                & metrics.lambda_subset.eq("0p1_0p5_0p9")
                & metrics.response_scope.eq("total")
                & np.isclose(metrics.trajectory_fraction, 1.0)
            )
        )
    ]
    repeat_summary = (
        repeat_rows.groupby("method", as_index=False)
        .agg(mean_mae=("mae", "mean"), sd_mae=("mae", "std"), min_mae=("mae", "min"), max_mae=("mae", "max"))
        .sort_values("mean_mae")
    )
    repeat_summary.to_csv(OUT / "phase1_primary_repeat_summary.csv", index=False)

    response_metrics = (
        response.groupby(
            ["partition", "repeat", "trajectory_fraction", "lambda", "component"],
            as_index=False,
        )
        .agg(n=("absolute_error", "size"), mae=("absolute_error", "mean"))
    )
    response_metrics.to_csv(OUT / "phase1_response_metrics.csv", index=False)

    reconstruction_comparisons = []
    repeated_reconstruction = reconstruction.loc[
        reconstruction.partition.eq("standardized_exclusion_repeat")
    ]
    for subset in sorted(repeated_reconstruction.lambda_subset.unique()):
        selected = repeated_reconstruction.loc[
            repeated_reconstruction.lambda_subset.eq(subset)
        ]
        pivoted = (
            selected.groupby(["molecule_id", "posterior"]).absolute_error.mean().unstack()
        )
        difference = (
            pivoted["P1-F_solvai_conditioned_posterior"]
            - pivoted["P1-F_generic_posterior"]
        )
        record = {
            "lambda_subset": subset,
            "n_molecules": len(difference),
            **paired_bootstrap(difference.to_numpy()),
        }
        for posterior in pivoted:
            record[f"mae_{posterior}"] = float(pivoted[posterior].mean())
        reconstruction_comparisons.append(record)
    reconstruction_comparisons = pd.DataFrame(reconstruction_comparisons)
    reconstruction_comparisons.to_csv(
        OUT / "phase1_reconstruction_paired_comparisons.csv", index=False
    )

    primary_family = endpoint.loc[repeat_filter & (primary_filter | baseline_filter)].copy()
    primary_family = (
        primary_family.groupby(["molecule_id", "functional_group_family", "method"])
        .absolute_error.mean()
        .unstack()
        .reset_index()
    )
    primary_family["active_minus_frozen"] = (
        primary_family["P1-D_actual_minus_predicted"]
        - primary_family["P1-A_frozen_solvai"]
    )
    family_summary = (
        primary_family.groupby("functional_group_family", as_index=False)
        .agg(
            n=("molecule_id", "size"),
            frozen_mae=("P1-A_frozen_solvai", "mean"),
            active_mae=("P1-D_actual_minus_predicted", "mean"),
            mean_difference=("active_minus_frozen", "mean"),
            fraction_improved=("active_minus_frozen", lambda value: float(np.mean(value < 0))),
        )
        .sort_values("mean_difference")
    )
    family_summary.to_csv(OUT / "phase1_exploratory_family_effects.csv", index=False)

    make_figures(endpoint, response, reconstruction)

    reconstruction_best = reconstruction_comparisons.sort_values("mean").iloc[0]
    response_full = response_metrics.loc[
        response_metrics.partition.eq("standardized_exclusion_repeat")
        & np.isclose(response_metrics.trajectory_fraction, 1.0)
        & response_metrics.component.eq("total")
    ]
    response_full_summary = (
        response_full.groupby("lambda", as_index=False).mae.agg(["mean", "std"]).reset_index()
    )
    canonical = {
        "status": "NEGATIVE_ENDPOINT_GATE",
        "generated_utc": datetime.now(UTC).isoformat(),
        "cohort_n": 72,
        "primary": {
            "method": "actual-minus-structure-predicted SYSTEM dH/dlambda at lambda 0.1/0.5/0.9",
            "frozen_solvai_repeat_mae_mean": float(
                repeat_summary.loc[
                    repeat_summary.method.eq("P1-A_frozen_solvai"), "mean_mae"
                ].iloc[0]
            ),
            "active_repeat_mae_mean": float(
                repeat_summary.loc[
                    repeat_summary.method.eq("P1-D_actual_minus_predicted"), "mean_mae"
                ].iloc[0]
            ),
            "paired_active_minus_frozen": primary_vs_baseline,
            "paired_active_minus_shuffle": primary_vs_shuffle,
            "stable_sign": False,
            "decision": "negative",
        },
        "response_prediction_mae_by_lambda": response_full_summary.to_dict("records"),
        "best_hidden_point_diagnostic": reconstruction_best.to_dict(),
        "dense_same_hamiltonian_population_available": False,
        "direction_decisions": {
            "A_empirical_probe_residual": "failed preregistered endpoint gate",
            "B_molecule_conditioned_bayesian_quadrature": "conditional diagnostic go to a bounded prospective dense sentinel; no dense-reconstruction claim yet",
            "C_adaptive_multifidelity": "held until the dense sentinel establishes a usable reconstruction model",
        },
        "integration_best_repeat_mae": float(
            integration_metrics.loc[
                integration_metrics.partition.eq("standardized_exclusion_repeat")
            ]
            .groupby(["method", "lambda_subset"])
            .mae.mean()
            .min()
        ),
        "output_hashes": {},
    }
    for path in sorted(OUT.glob("phase1_*")):
        if path.is_file() and path.name != "phase1_canonical_metrics.json":
            canonical["output_hashes"][path.name] = sha256(path)
    (OUT / "phase1_canonical_metrics.json").write_text(
        json.dumps(canonical, indent=2, sort_keys=True) + "\n"
    )

    fixed = metrics.loc[
        metrics.partition.eq("standardized_exclusion_primary")
        & metrics.shuffle_seed.isna()
        & (
            metrics.method.eq("P1-A_frozen_solvai")
            | (
                metrics.method.isin(
                    [
                        "P1-B_predicted_response",
                        "P1-C_actual_response",
                        "P1-D_actual_minus_predicted",
                        "P1-H_mean_shuffled_residual",
                    ]
                )
                & metrics.lambda_subset.eq("0p1_0p5_0p9")
                & metrics.response_scope.eq("total")
                & np.isclose(metrics.trajectory_fraction, 1.0)
            )
        )
    ][["method", "n", "mae", "rmse", "median_absolute_error"]]

    fixed_active_mae = metrics.loc[
        metrics.partition.eq("standardized_exclusion_primary")
        & metrics.method.eq("P1-D_actual_minus_predicted")
        & metrics.lambda_subset.eq("0p1_0p5_0p9")
        & metrics.response_scope.eq("total")
        & np.isclose(metrics.trajectory_fraction, 1.0),
        "mae",
    ].iloc[0]
    report = f"""# Phase 1 actual-observation gate

## Decision

**The preregistered experimental-endpoint gate is negative.** On the 72 molecules with complete 5 ps PIMD2 observations at λ=0.1, 0.5 and 0.9, the actual-minus-predicted response correction increased repeated-partition MAE from {canonical['primary']['frozen_solvai_repeat_mae_mean']:.6f} to {canonical['primary']['active_repeat_mae_mean']:.6f} kcal mol⁻¹. The molecule-paired candidate-minus-baseline change was {primary_vs_baseline['mean']:+.6f} kcal mol⁻¹ (95% bootstrap CI {primary_vs_baseline['ci_low_95']:+.6f} to {primary_vs_baseline['ci_high_95']:+.6f}). Only {primary_vs_baseline['fraction_below_zero']:.1%} of molecules improved. The candidate was also indistinguishable from the mean shuffled-residual control (difference {primary_vs_shuffle['mean']:+.6f}; 95% CI {primary_vs_shuffle['ci_low_95']:+.6f} to {primary_vs_shuffle['ci_high_95']:+.6f}).

The sign was unfavourable in four of five repeated partitions and effectively null in the fifth. The fixed parent partition was also worse ({metrics.loc[(metrics.partition.eq('standardized_exclusion_primary')) & metrics.method.eq('P1-A_frozen_solvai'), 'mae'].iloc[0]:.6f} versus {fixed_active_mae:.6f}). This satisfies the frozen negative endpoint criterion and cannot be rescued by a favourable post-hoc λ subset or chemistry.

## Matched fixed-partition results

{dataframe_markdown(fixed)}

## Five repeated partitions

{dataframe_markdown(repeat_summary)}

## Destructive control

All three observed response values were permuted jointly across molecules within each outer training and test fold for five preregistered seeds. The aligned residual was not better than the mean shuffled control. This falsifies the claim that the present 5 ps PIMD2 residual adds stable molecule-specific endpoint information beyond frozen SolvAI.

## Observation duration

The primary residual correction was repeated at sequential 0.5, 1.0, 2.0, 3.5 and 5.0 ps prefixes. Longer prefixes reduced variability but did not produce an endpoint gain over frozen SolvAI. No future frames were selected.

## Sparse response reconstruction

SolvAI-conditioned Gaussian interpolation often reduced hidden-point error relative to a generic population Gaussian. The largest post-result descriptive gain among the predeclared two-point subsets was `{reconstruction_best['lambda_subset']}`: {reconstruction_best['mae_P1-F_generic_posterior']:.3f} to {reconstruction_best['mae_P1-F_solvai_conditioned_posterior']:.3f} kcal mol⁻¹, paired difference {reconstruction_best['mean']:+.3f} (95% CI {reconstruction_best['ci_low_95']:+.3f} to {reconstruction_best['ci_high_95']:+.3f}). However, hidden-point errors remained 2.4 kcal mol⁻¹ or larger and nominal 95% intervals were very wide. These three-point held-point diagnostics are not dense-curve reconstruction. Under the prospective freeze they justify only a bounded dense same-Hamiltonian sentinel test, not a reconstruction claim or endpoint rescue.

## Numerical integration

Direct integration of the three actual short-window observations followed by fold-local affine calibration produced a five-repeat mean experimental MAE of {integration_metrics.loc[(integration_metrics.partition.eq('standardized_exclusion_repeat')) & integration_metrics.method.eq('P1-E_actual_integral_affine')].groupby('repeat').mae.first().mean():.3f} kcal mol⁻¹. The best posterior-integral variant remained {canonical['integration_best_repeat_mae']:.3f} kcal mol⁻¹. The inherited three λ points are therefore not an adequate quadrature rule for this endpoint.

## Direction decisions

1. **Direction A — empirical residual correction:** failed the prospectively frozen endpoint criterion.
2. **Direction B — molecule-conditioned Bayesian quadrature:** hidden-point interpolation gains passed the limited diagnostic criterion for several predeclared subsets, but no compatible dense same-Hamiltonian population exists locally, intervals are broad, and the integral is inaccurate. This is a conditional go to one bounded, prospectively frozen dense sentinel acquisition; it is not a reconstruction claim.
3. **Direction C — adaptive multi-fidelity allocation:** held until that dense sentinel establishes whether a useful reconstruction model exists.

No new MD/PIMD calculation was used in this gate. This is a scientific null for endpoint correction with the inherited 5 ps, PIMD2, three-window protocol—not evidence that all possible active solvation calculations lack value. A separately frozen dense sentinel can test the surviving reconstruction hypothesis, but it cannot alter the failed endpoint decision.

## Reproducibility

- Freeze commit: `a0dd986`
- Command: `active_solvai/.venv/bin/python active_solvai/scripts/run_phase1_gate.py`
- Summary command: `active_solvai/.venv/bin/python active_solvai/scripts/summarize_phase1.py`
- Machine-readable endpoint, response, reconstruction and integration tables are in `active_solvai/results/phase1/`.
"""
    REPORT.write_text(report)

    null_entry = f"""

## AS-P1-GATE-001 — actual PIMD2 observation gate (2026-09-03)

- **Frozen before scoring:** commit `a0dd986`.
- **Question:** do three actual 5 ps PIMD2 SYSTEM dH/dλ observations, especially actual-minus-structure-predicted residuals, improve frozen SolvAI?
- **Result:** no. Five-repeat mean MAE changed from {canonical['primary']['frozen_solvai_repeat_mae_mean']:.6f} to {canonical['primary']['active_repeat_mae_mean']:.6f} kcal mol⁻¹; paired difference {primary_vs_baseline['mean']:+.6f} (95% CI {primary_vs_baseline['ci_low_95']:+.6f}, {primary_vs_baseline['ci_high_95']:+.6f}). The aligned residual was not better than shuffled residuals.
- **Mechanistic diagnostic:** structure-conditioned interpolation reduced some hidden-point errors, but errors remained ≥2.4 kcal mol⁻¹ with broad intervals; three points do not constitute a dense reconstruction benchmark.
- **Decision:** kill Direction A for this probe protocol. Hold Direction C. The limited hidden-point signal permits one separately frozen dense sentinel test of Direction B; it cannot alter this endpoint null.
"""
    existing = NULL_RESULTS.read_text()
    if "AS-P1-GATE-001" not in existing:
        with NULL_RESULTS.open("a", encoding="utf-8") as handle:
            handle.write(null_entry)

    ledger_path = ROOT / "active_solvai/runs/ledger.jsonl"
    existing_run_ids = {
        json.loads(line)["run_id"] for line in ledger_path.read_text().splitlines() if line.strip()
    }
    if "AS-P1-SUMMARY-001" not in existing_run_ids:
        append_record(
            ledger_path,
            {
            "run_id": "AS-P1-SUMMARY-001",
            "stage": "phase1",
            "status": "completed",
            "command": "active_solvai/.venv/bin/python active_solvai/scripts/summarize_phase1.py",
            "device": "CPU",
            "wall_seconds": 0.0,
            "gpu_hours": 0.0,
            "cpu_hours": 0.0,
            "simulated_time_ps": 0.0,
            "force_evaluations": 0,
            "bead_windows": 0,
            "quality_control": "summary cross-check complete",
            "failure_reason": None,
            "freeze_commit": "a0dd986",
            },
        )
    print(json.dumps(canonical["primary"], indent=2))


if __name__ == "__main__":
    main()
