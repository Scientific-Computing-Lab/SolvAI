"""Summarize the frozen prospective dense-sentinel replay."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from active_solvai.ledger import append_record, sha256

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"
RESULTS = ACTIVE_ROOT / "results/phase2"
FIGURES = ACTIVE_ROOT / "figures/phase2"
REPORT = ACTIVE_ROOT / "reports/DENSE_SENTINEL_REPLAY.md"

DEPLOYABLE_COMPARATORS = (
    "fixed_solvai_bq",
    "uniform_solvai_bq",
    "random_solvai_bq",
    "curvature_solvai_bq",
    "generic_bq",
    "fixed_direct",
    "uniform_direct",
)


def fmt(value: float) -> str:
    return "—" if not np.isfinite(value) else f"{value:.3f}"


def main() -> None:
    metrics = pd.read_csv(RESULTS / "dense_replay_metrics.csv")
    comparisons = pd.read_csv(RESULTS / "dense_replay_paired_comparisons.csv")
    predictions = pd.read_parquet(RESULTS / "dense_replay_predictions.parquet")
    FIGURES.mkdir(parents=True, exist_ok=True)

    primary_rows: list[dict[str, object]] = []
    condition_one = False
    for budget in (5, 7):
        candidates = metrics.loc[
            metrics.total_windows.eq(budget) & metrics.method.isin(DEPLOYABLE_COMPARATORS)
        ]
        best = candidates.sort_values("integral_mae_kcal_mol").iloc[0]
        comparison = comparisons.loc[
            comparisons.total_windows.eq(budget) & comparisons.comparator.eq(best.method)
        ].iloc[0]
        active = metrics.loc[
            metrics.total_windows.eq(budget) & metrics.method.eq("active_solvai_bq")
        ].iloc[0]
        passed = bool(
            comparison.candidate_minus_comparator_mae <= -0.10
            and comparison.ci90_high < 0.0
            and comparison.fraction_candidate_improved >= 0.75
            and 0.75 <= active.coverage_90 <= 1.0
        )
        condition_one |= passed
        primary_rows.append(
            {
                "budget": budget,
                "active_mae": float(active.integral_mae_kcal_mol),
                "active_coverage_90": float(active.coverage_90),
                "strongest_comparator": str(best.method),
                "comparator_mae": float(best.integral_mae_kcal_mol),
                "paired_difference": float(comparison.candidate_minus_comparator_mae),
                "ci90_low": float(comparison.ci90_low),
                "ci90_high": float(comparison.ci90_high),
                "fraction_improved": float(comparison.fraction_candidate_improved),
                "condition_one_passed": passed,
            }
        )

    candidate_tolerance = metrics.loc[
        metrics.method.eq("active_solvai_bq") & (metrics.integral_mae_kcal_mol <= 0.20)
    ]
    candidate_min = (
        int(candidate_tolerance.total_windows.min()) if len(candidate_tolerance) else None
    )
    comparator_mins: dict[str, int | None] = {}
    for method in DEPLOYABLE_COMPARATORS:
        rows = metrics.loc[metrics.method.eq(method) & (metrics.integral_mae_kcal_mol <= 0.20)]
        comparator_mins[method] = int(rows.total_windows.min()) if len(rows) else None
    finite_comparator_mins = [value for value in comparator_mins.values() if value is not None]
    condition_two = bool(
        candidate_min is not None
        and finite_comparator_mins
        and candidate_min <= 0.8 * min(finite_comparator_mins)
    )
    decision = "positive" if condition_one or condition_two else "negative"

    plot_methods = [
        "active_solvai_bq",
        "generic_bq",
        "fixed_solvai_bq",
        "uniform_solvai_bq",
        "random_solvai_bq",
        "fixed_direct",
        "oracle_non_deployable",
    ]
    colors = {
        "active_solvai_bq": "#0072B2",
        "generic_bq": "#999999",
        "fixed_solvai_bq": "#E69F00",
        "uniform_solvai_bq": "#56B4E9",
        "random_solvai_bq": "#CC79A7",
        "fixed_direct": "#009E73",
        "oracle_non_deployable": "#D55E00",
    }
    labels = {method: method.replace("_", " ") for method in plot_methods}
    fig, axis = plt.subplots(figsize=(7.2, 4.5))
    for method in plot_methods:
        rows = metrics.loc[metrics.method.eq(method)].sort_values("total_windows")
        axis.plot(
            rows.total_windows,
            rows.integral_mae_kcal_mol,
            marker="o",
            linewidth=1.8,
            markersize=4,
            color=colors[method],
            label=labels[method],
        )
    axis.set_xlabel("Observed PIMD2 windows (including three inherited windows)")
    axis.set_ylabel("Integral MAE to dense PIMD2 (kcal mol$^{-1}$)")
    axis.legend(frameon=False, fontsize=8, ncol=2)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(FIGURES / f"dense_cost_error_frontier.{suffix}", dpi=300)
    plt.close(fig)

    selected = predictions.loc[
        predictions.method.isin(["active_solvai_bq", "generic_bq", "fixed_solvai_bq"])
        & predictions.total_windows.eq(5)
        & predictions.schedule_replicate.eq(0)
    ]
    pivot = selected.pivot(
        index="molecule_name", columns="method", values="absolute_integral_error_kcal_mol"
    ).sort_values("active_solvai_bq")
    fig, axis = plt.subplots(figsize=(7.2, 4.5))
    locations = np.arange(len(pivot))
    width = 0.24
    for offset, method in zip((-width, 0.0, width), pivot.columns, strict=True):
        axis.bar(
            locations + offset, pivot[method], width, color=colors[method], label=labels[method]
        )
    axis.set_xticks(locations, pivot.index, rotation=40, ha="right")
    axis.set_ylabel("Absolute integral error (kcal mol$^{-1}$)")
    axis.legend(frameon=False, fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(FIGURES / f"dense_molecule_errors_budget5.{suffix}", dpi=300)
    plt.close(fig)

    table = metrics.loc[
        metrics.total_windows.isin([3, 5, 7, 9, 15]) & metrics.method.isin(plot_methods)
    ].copy()
    lines = [
        "# Prospective dense PIMD2 sentinel replay",
        "",
        f"**Frozen decision: {decision.upper()}.**",
        "",
        "This is a prospective test of short-window, same-Hamiltonian PIMD2 response reconstruction. It is not a test of experimental endpoint improvement, full PIMD8 convergence or blind chemistry.",
        "",
        "## Primary frozen comparisons",
        "",
        "| Total windows | Active MAE | 90% coverage | Strongest comparator | Comparator MAE | Paired difference | 90% CI | Improved | Pass |",
        "|---:|---:|---:|---|---:|---:|---:|---:|:---:|",
    ]
    for row in primary_rows:
        lines.append(
            f"| {row['budget']} | {row['active_mae']:.3f} | {row['active_coverage_90']:.3f} | "
            f"{row['strongest_comparator']} | {row['comparator_mae']:.3f} | "
            f"{row['paired_difference']:+.3f} | [{row['ci90_low']:+.3f}, {row['ci90_high']:+.3f}] | "
            f"{row['fraction_improved']:.3f} | {'yes' if row['condition_one_passed'] else 'no'} |"
        )
    lines.extend(
        [
            "",
            "## Accuracy–cost frontier",
            "",
            "| Method | Windows | Integral MAE | Hidden-curve MAE | 90% coverage | Mean 90% width | Mean measured wall time (s) |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in table.sort_values(["total_windows", "integral_mae_kcal_mol"]).to_dict("records"):
        lines.append(
            f"| {row['method']} | {int(row['total_windows'])} | {row['integral_mae_kcal_mol']:.3f} | "
            f"{row['hidden_curve_mae_kcal_mol']:.3f} | {fmt(row['coverage_90'])} | "
            f"{fmt(row['mean_interval_width_90_kcal_mol'])} | {row['mean_measured_window_wall_seconds']:.1f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "The frozen positive criterion was met; Direction B survives to a separately frozen multi-fidelity test."
                if decision == "positive"
                else "The frozen positive criterion was not met. Direction B is killed for this protocol and Direction C is not launched."
            ),
            "",
            f"The active policy first reached 0.20 kcal mol⁻¹ MAE at {candidate_min if candidate_min is not None else 'no tested'} windows. Comparator minima were `{json.dumps(comparator_mins, sort_keys=True)}`.",
            "",
            "The oracle is non-deployable and cannot support a method claim. All twelve molecules were already development-exposed in the parent project; only their twelve added dense-window responses were prospective.",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n")
    canonical = {
        "schema_version": 1,
        "decision": decision,
        "condition_one_passed": condition_one,
        "condition_two_passed": condition_two,
        "primary_comparisons": primary_rows,
        "active_minimum_windows_below_0p20": candidate_min,
        "comparator_minimum_windows_below_0p20": comparator_mins,
        "inputs": {
            "predictions_sha256": sha256(RESULTS / "dense_replay_predictions.parquet"),
            "metrics_sha256": sha256(RESULTS / "dense_replay_metrics.csv"),
            "comparisons_sha256": sha256(RESULTS / "dense_replay_paired_comparisons.csv"),
        },
    }
    canonical_path = RESULTS / "dense_replay_canonical_metrics.json"
    canonical_path.write_text(json.dumps(canonical, indent=2) + "\n")
    ledger = ACTIVE_ROOT / "runs/ledger.jsonl"
    if not any(
        json.loads(line).get("run_id") == "AS-P2-REPLAY-001"
        for line in ledger.read_text().splitlines()
        if line.strip()
    ):
        append_record(
            ledger,
            {
                "run_id": "AS-P2-REPLAY-001",
                "timestamp_utc": datetime.now(UTC).isoformat(),
                "stage": "phase2_dense_prospective",
                "action": "frozen dense common-pool replay and summary",
                "status": "completed",
                "quality_control": decision,
                "device": "CPU",
                "gpu_hours": 0.0,
                "bead_windows": 0,
                "simulated_time_ps": 0.0,
                "force_evaluations": 0,
                "output": str(canonical_path),
                "output_sha256": sha256(canonical_path),
            },
        )
    print(json.dumps(canonical, indent=2))


if __name__ == "__main__":
    main()
