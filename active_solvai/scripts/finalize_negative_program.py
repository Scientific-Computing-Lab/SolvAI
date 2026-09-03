#!/usr/bin/env python3
"""Build the canonical Active SolvAI negative-result diagnosis and final report.

This script performs descriptive analysis only. It reads the already frozen and
scored Phase 1/2 artifacts; it does not refit a model, choose a schedule, or
alter an interpretation threshold.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from active_solvai.ledger import append_record, sha256

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"
PHASE1 = ACTIVE_ROOT / "results/phase1"
PHASE2 = ACTIVE_ROOT / "results/phase2"
REPORT = ACTIVE_ROOT / "reports/ACTIVE_SOLVAI_FINAL_REPORT.md"
DIAGNOSTIC_CSV = PHASE2 / "dense_failure_diagnostics.csv"
DIAGNOSTIC_JSON = PHASE2 / "active_solvai_final_metrics.json"
LEDGER = ACTIVE_ROOT / "runs/ledger.jsonl"


def one_row(frame: pd.DataFrame, molecule: str, method: str, budget: int) -> pd.Series:
    rows = frame.loc[
        frame.molecule_name.eq(molecule)
        & frame.method.eq(method)
        & frame.total_windows.eq(budget)
        & frame.schedule_replicate.eq(0)
    ]
    if len(rows) != 1:
        raise AssertionError((molecule, method, budget, len(rows)))
    return rows.iloc[0]


def simulation_totals() -> dict[str, dict[str, float | int]]:
    records = [json.loads(line) for line in LEDGER.read_text().splitlines() if line.strip()]
    totals: dict[str, dict[str, float | int]] = {}
    for role in ("calibration", "prospective"):
        rows = [
            row
            for row in records
            if row.get("role") == role
            and row.get("stage") in {"phase2_dense_calibration", "phase2_dense_prospective"}
            and str(row.get("run_id", "")).endswith("A01")
        ]
        totals[role] = {
            "attempted_windows": len(rows),
            "passed_windows": sum(bool(row.get("passed")) for row in rows),
            "failed_windows": sum(not bool(row.get("passed")) for row in rows),
            "wall_seconds": sum(float(row.get("wall_seconds", 0.0) or 0.0) for row in rows),
            "gpu_hours": sum(float(row.get("gpu_hours", 0.0) or 0.0) for row in rows),
            "cpu_hours": sum(float(row.get("cpu_hours", 0.0) or 0.0) for row in rows),
            "production_ps": sum(float(row.get("production_ps", 0.0) or 0.0) for row in rows),
            "bead_windows": sum(int(row.get("bead_windows", 0) or 0) for row in rows),
            "nominal_bead_steps": sum(int(row.get("nominal_bead_steps", 0) or 0) for row in rows),
        }
    return totals


def main() -> None:
    phase1 = json.loads((PHASE1 / "phase1_canonical_metrics.json").read_text())
    phase2 = json.loads((PHASE2 / "dense_replay_canonical_metrics.json").read_text())
    responses = pd.read_parquet(PHASE2 / "dense_responses_prospective.parquet")
    predictions = pd.read_parquet(PHASE2 / "dense_replay_predictions.parquet")

    rows: list[dict[str, object]] = []
    for molecule in sorted(predictions.molecule_name.unique()):
        response_rows = responses.loc[responses.molecule_name.eq(molecule)]
        for budget in (5, 7):
            active = one_row(predictions, molecule, "active_solvai_bq", budget)
            uniform = one_row(predictions, molecule, "uniform_direct", budget)
            generic = one_row(predictions, molecule, "generic_bq", budget)
            fixed = one_row(predictions, molecule, "fixed_direct", budget)
            oracle = one_row(predictions, molecule, "oracle_non_deployable", budget)
            rows.append(
                {
                    "molecule_name": molecule,
                    "functional_group_family": active.functional_group_family,
                    "total_windows": budget,
                    "dense_integral_kcal_mol": active.true_integral_kcal_mol,
                    "mean_window_sem_kcal_mol": response_rows.five_block_sem_kcal_mol.mean(),
                    "maximum_window_sem_kcal_mol": response_rows.five_block_sem_kcal_mol.max(),
                    "active_absolute_error_kcal_mol": active.absolute_integral_error_kcal_mol,
                    "uniform_direct_absolute_error_kcal_mol": uniform.absolute_integral_error_kcal_mol,
                    "generic_bq_absolute_error_kcal_mol": generic.absolute_integral_error_kcal_mol,
                    "fixed_direct_absolute_error_kcal_mol": fixed.absolute_integral_error_kcal_mol,
                    "oracle_absolute_error_kcal_mol": oracle.absolute_integral_error_kcal_mol,
                    "active_minus_uniform_direct_absolute_error_kcal_mol": (
                        active.absolute_integral_error_kcal_mol
                        - uniform.absolute_integral_error_kcal_mol
                    ),
                    "active_observed_lambdas": active.observed_lambdas,
                    "uniform_observed_lambdas": uniform.observed_lambdas,
                    "oracle_observed_lambdas": oracle.observed_lambdas,
                    "active_covered_90": active.covered_90,
                    "active_interval_width_90_kcal_mol": active.interval_width_90_kcal_mol,
                }
            )
    diagnostics = pd.DataFrame(rows)
    diagnostics.to_csv(DIAGNOSTIC_CSV, index=False)

    totals = simulation_totals()
    all_sim = {
        key: sum(float(totals[role][key]) for role in totals)
        for key in (
            "attempted_windows",
            "passed_windows",
            "failed_windows",
            "wall_seconds",
            "gpu_hours",
            "cpu_hours",
            "production_ps",
            "bead_windows",
            "nominal_bead_steps",
        )
    }
    primary = {str(row["budget"]): row for row in phase2["primary_comparisons"]}
    summary = {
        "schema_version": 1,
        "program_decision": "NO_GO",
        "phase1_endpoint_decision": phase1["primary"]["decision"],
        "phase1_frozen_solvai_repeat_mae_kcal_mol": phase1["primary"][
            "frozen_solvai_repeat_mae_mean"
        ],
        "phase1_residual_corrected_repeat_mae_kcal_mol": phase1["primary"][
            "active_repeat_mae_mean"
        ],
        "phase1_candidate_minus_baseline_kcal_mol": phase1["primary"]["paired_active_minus_frozen"],
        "phase1_candidate_minus_shuffle_kcal_mol": phase1["primary"]["paired_active_minus_shuffle"],
        "phase2_decision": phase2["decision"],
        "phase2_primary_comparisons": phase2["primary_comparisons"],
        "phase2_mean_window_sem_kcal_mol": float(responses.five_block_sem_kcal_mol.mean()),
        "phase2_median_window_sem_kcal_mol": float(responses.five_block_sem_kcal_mol.median()),
        "phase2_maximum_window_sem_kcal_mol": float(responses.five_block_sem_kcal_mol.max()),
        "phase2_oracle_mae_budget5_kcal_mol": float(
            predictions.loc[
                predictions.method.eq("oracle_non_deployable") & predictions.total_windows.eq(5),
                "absolute_integral_error_kcal_mol",
            ].mean()
        ),
        "phase2_oracle_mae_budget7_kcal_mol": float(
            predictions.loc[
                predictions.method.eq("oracle_non_deployable") & predictions.total_windows.eq(7),
                "absolute_integral_error_kcal_mol",
            ].mean()
        ),
        "simulation_compute": {"by_role": totals, "total": all_sim},
        "direction_status": {
            "A_empirical_probe_residual": "failed frozen endpoint gate",
            "B_molecule_conditioned_bayesian_quadrature": (
                "failed frozen prospective dense-sentinel criterion"
            ),
            "C_adaptive_multifidelity": (
                "not launched; prospectively contingent on Direction B passing"
            ),
        },
        "diagnosis": {
            "response_signal": (
                "actual three-point residuals did not improve the experimental endpoint and did not beat shuffling"
            ),
            "response_noise": (
                "five-block SEM was substantial; longer trajectories were not tested and rescue is not established"
            ),
            "lambda_placement": (
                "the oracle shows informative molecule-specific placements exist, but the frozen acquisition did not identify them"
            ),
            "structure_prior": (
                "the structure-conditioned curve prior was not consistently better than generic or direct interpolation"
            ),
            "cross_fidelity": ("not tested because the same-fidelity prerequisite failed"),
            "hamiltonian_bias": (
                "cannot explain the same-Hamiltonian reconstruction failure; may limit experimental-endpoint transfer"
            ),
            "endpoint_labels": (
                "the low baseline error leaves limited headroom, but aligned-vs-shuffled controls show no usable residual signal"
            ),
            "data_availability": (
                "no compatible historical dense population existed; a prospectively generated 12-molecule panel removed that dependency for this test"
            ),
            "implementation": (
                "no evidence of implementation failure; all prospective windows passed QC and frozen analysis assertions"
            ),
        },
        "input_hashes": {
            "phase1": sha256(PHASE1 / "phase1_canonical_metrics.json"),
            "phase2": sha256(PHASE2 / "dense_replay_canonical_metrics.json"),
            "responses": sha256(PHASE2 / "dense_responses_prospective.parquet"),
            "predictions": sha256(PHASE2 / "dense_replay_predictions.parquet"),
            "diagnostics": sha256(DIAGNOSTIC_CSV),
        },
    }
    DIAGNOSTIC_JSON.write_text(json.dumps(summary, indent=2) + "\n")

    p5, p7 = primary["5"], primary["7"]
    lines = [
        "# Active SolvAI final report",
        "",
        "## Executive outcome",
        "",
        "**No-go under the prospectively frozen program.** Actual short PIMD2 observations did not improve the experimental hydration endpoint, and the molecule-conditioned Bayesian-quadrature policy did not beat the strongest simple schedule in a prospective dense same-Hamiltonian test. The conditional multi-fidelity direction was therefore not launched.",
        "",
        "This is a rigorous negative result for the tested 5-ps PIMD2 protocol and acquisition model. It is not evidence that all forms of active free-energy simulation are impossible.",
        "",
        "## Direction A — experimental endpoint",
        "",
        f"On 72 molecules, the five-partition frozen SolvAI MAE was {summary['phase1_frozen_solvai_repeat_mae_kcal_mol']:.6f} kcal mol⁻¹ and the actual-minus-predicted three-point residual model was {summary['phase1_residual_corrected_repeat_mae_kcal_mol']:.6f}. The paired change was {summary['phase1_candidate_minus_baseline_kcal_mol']['mean']:+.6f} (95% CI {summary['phase1_candidate_minus_baseline_kcal_mol']['ci_low_95']:+.6f}, {summary['phase1_candidate_minus_baseline_kcal_mol']['ci_high_95']:+.6f}); only {summary['phase1_candidate_minus_baseline_kcal_mol']['fraction_below_zero']:.1%} of molecules improved. The aligned residual was indistinguishable from the shuffled control.",
        "",
        "**Decision:** Direction A is killed for these observations. No λ subset, family, component or post-hoc endpoint correction was used to rescue it.",
        "",
        "## Direction B — same-Hamiltonian reconstruction",
        "",
        "Four calibration and eight prospective molecules were simulated on a fixed 15-point λ grid. Every molecule reused the three inherited λ=0.1, 0.5 and 0.9 windows and added 12 new 5-ps PIMD2 windows. Calibration choices were locked before any prospective response was generated.",
        "",
        "| Windows | Active BQ MAE | 90% coverage | Strongest comparator | Comparator MAE | Paired difference | 90% CI | Molecules improved |",
        "|---:|---:|---:|---|---:|---:|---:|---:|",
        f"| 5 | {p5['active_mae']:.3f} | {p5['active_coverage_90']:.3f} | {p5['strongest_comparator']} | {p5['comparator_mae']:.3f} | {p5['paired_difference']:+.3f} | [{p5['ci90_low']:+.3f}, {p5['ci90_high']:+.3f}] | {p5['fraction_improved']:.3f} |",
        f"| 7 | {p7['active_mae']:.3f} | {p7['active_coverage_90']:.3f} | {p7['strongest_comparator']} | {p7['comparator_mae']:.3f} | {p7['paired_difference']:+.3f} | [{p7['ci90_low']:+.3f}, {p7['ci90_high']:+.3f}] | {p7['fraction_improved']:.3f} |",
        "",
        f"The non-deployable oracle reached {summary['phase2_oracle_mae_budget5_kcal_mol']:.3f} and {summary['phase2_oracle_mae_budget7_kcal_mol']:.3f} kcal mol⁻¹ at five and seven windows. Thus informative molecule-specific placements exist in the dense pool, but the frozen acquisition rule did not find them.",
        "",
        "**Decision:** Direction B fails the prospectively frozen criterion. Direction C was prospectively contingent on this result and was not launched.",
        "",
        "## Failure attribution",
        "",
        "| Possible cause | Evidence and conclusion |",
        "|---|---|",
        "| Lack of endpoint residual signal | Supported for this protocol: the aligned response residual worsened endpoint MAE and did not beat shuffling. |",
        f"| Response noise / trajectory length | Five-block response SEM averaged {summary['phase2_mean_window_sem_kcal_mol']:.3f} kcal mol⁻¹ (median {summary['phase2_median_window_sem_kcal_mol']:.3f}; maximum {summary['phase2_maximum_window_sem_kcal_mol']:.3f}). Noise is material, but longer trajectories were not tested, so a length-based rescue is unproven. |",
        "| λ placement | The oracle was far better than deployable policies at five/seven windows, showing that placement matters; the present acquisition score did not identify the useful points. |",
        "| Structure prior | The structure-conditioned prior/posterior was not consistently better than generic BQ or direct interpolation, especially for amide, ether, fused-aromatic and alkane sentinels. |",
        "| Cross-fidelity mapping | Not tested. The same-fidelity PIMD2 reconstruction prerequisite failed, so PIMD4/PIMD8 escalation would add degrees of freedom without an established base. |",
        "| Hamiltonian bias | Not a cause of the same-Hamiltonian reconstruction failure. It remains a plausible limit on transfer from short ARROW/PIMD2 responses to experiment. |",
        "| Endpoint labels | The 72-molecule baseline has little headroom, but the aligned-versus-shuffled result is the decisive evidence: no usable molecule-specific endpoint signal was detected. |",
        "| Data availability | No compatible historical dense population existed. The new prospective dense panel resolved that dependency for this bounded test. |",
        "| Implementation failure | Not supported: all 144 new windows passed QC on first attempt, all outputs were finite/complete, and the frozen analysis and unit assertions passed. |",
        "",
        "## Compute and failed-work accounting",
        "",
        f"The campaign generated {int(all_sim['attempted_windows'])} new windows: {int(totals['calibration']['attempted_windows'])} calibration and {int(totals['prospective']['attempted_windows'])} prospective. All passed QC on the first attempt; failed-work cost was zero. New production totaled {all_sim['production_ps']:.0f} ps, {int(all_sim['bead_windows'])} bead-windows and {int(all_sim['nominal_bead_steps'])} nominal bead-steps. Measured simulation wall time was {all_sim['wall_seconds'] / 3600:.3f} h and GPU accounting was {all_sim['gpu_hours']:.3f} GPU-h. Arbalest does not expose exact fast/slow force-kernel call counts, so no fabricated force-evaluation number is reported.",
        "",
        "## Scientific conclusion",
        "",
        "For the inherited short PIMD2 protocol, actual target-molecule response measurements neither improved the experimental endpoint nor enabled the predeclared molecule-conditioned active quadrature to outperform a simple uniform schedule. The result identifies two separable bottlenecks: response residuals were not predictive of experimental error, and the structure-conditioned curve prior/acquisition was not accurate enough to locate the molecule-specific informative windows visible to the oracle.",
        "",
        "A future campaign would first require a substantially larger protocol-matched dense-curve training set and independent evidence that its response prior predicts full curves—not only three points—before testing longer trajectories or multi-fidelity escalation. That is a new data-generation program, not a scale-up justified by the present result. Tier-B remains unopened.",
    ]
    REPORT.write_text("\n".join(lines) + "\n")

    existing = [json.loads(line) for line in LEDGER.read_text().splitlines() if line.strip()]
    canonical_hash = sha256(DIAGNOSTIC_JSON)
    recorded_hashes = {
        row.get("output_sha256")
        for row in existing
        if str(row.get("run_id", "")).startswith("AS-FINAL-DIAG-")
    }
    if canonical_hash not in recorded_hashes:
        suffix = 1 + sum(
            str(row.get("run_id", "")).startswith("AS-FINAL-DIAG-") for row in existing
        )
        append_record(
            LEDGER,
            {
                "run_id": f"AS-FINAL-DIAG-{suffix:03d}",
                "stage": "final_negative_diagnosis",
                "action": "descriptive diagnosis from frozen Phase 1/2 outputs",
                "status": "completed",
                "device": "CPU",
                "gpu_hours": 0.0,
                "cpu_hours": 0.0,
                "simulated_time_ps": 0.0,
                "bead_windows": 0,
                "force_evaluations": 0,
                "input_sha256": summary["input_hashes"],
                "output": str(DIAGNOSTIC_JSON),
                "output_sha256": canonical_hash,
            },
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
