#!/usr/bin/env python3
"""Freeze the Phase 1 confirmatory results into tables and a scientific report."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "confirmatory"
REPORT = ROOT / "reports" / "CONFIRMATORY_ANALYSIS.md"
BOOTSTRAP_SEED = 20260828
BOOTSTRAP_REPLICATES = 100_000


def paired_bootstrap(difference: np.ndarray) -> dict[str, float]:
    """Bootstrap a molecule-level candidate-minus-reference error difference."""
    difference = np.asarray(difference, dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n = len(difference)
    chunk = 2_000
    means: list[np.ndarray] = []
    for start in range(0, BOOTSTRAP_REPLICATES, chunk):
        size = min(chunk, BOOTSTRAP_REPLICATES - start)
        indices = rng.integers(0, n, size=(size, n))
        means.append(difference[indices].mean(axis=1))
    samples = np.concatenate(means)
    return {
        "difference": float(difference.mean()),
        "ci_low_95": float(np.quantile(samples, 0.025)),
        "ci_high_95": float(np.quantile(samples, 0.975)),
        "probability_lower_error": float(np.mean(samples < 0.0)),
        "fraction_molecules_improved": float(np.mean(difference < 0.0)),
    }


def metrics(truth: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    residual = np.asarray(truth, dtype=float) - np.asarray(prediction, dtype=float)
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((np.asarray(truth) - np.mean(truth)) ** 2))
    return {
        "n": len(residual),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "median_absolute_error": float(np.median(np.abs(residual))),
        "r2": float(1.0 - ss_res / ss_tot),
    }


def prediction_vector(frame: pd.DataFrame, method: str) -> pd.DataFrame:
    selected = frame.loc[frame.method.eq(method)].copy()
    return selected.sort_values("molecule_id").reset_index(drop=True)


def comparison_row(
    analysis: str,
    candidate: pd.DataFrame,
    reference: pd.DataFrame,
    candidate_label: str,
    reference_label: str,
) -> dict[str, object]:
    left = prediction_vector(candidate, candidate_label).set_index("molecule_id")
    right = prediction_vector(reference, reference_label).set_index("molecule_id")
    if set(left.index) != set(right.index):
        raise AssertionError(f"Molecule mismatch in {analysis}")
    left = left.loc[sorted(left.index)]
    right = right.loc[left.index]
    if not np.allclose(left.y_true, right.y_true, rtol=0.0, atol=1e-12):
        raise AssertionError(f"Truth mismatch in {analysis}")
    candidate_error = np.abs(left.y_true.to_numpy() - left.y_pred.to_numpy())
    reference_error = np.abs(right.y_true.to_numpy() - right.y_pred.to_numpy())
    result: dict[str, object] = {
        "analysis": analysis,
        "candidate": candidate_label,
        "reference": reference_label,
        "candidate_mae": float(candidate_error.mean()),
        "reference_mae": float(reference_error.mean()),
    }
    result.update(paired_bootstrap(candidate_error - reference_error))
    result["outcome"] = (
        "positive"
        if result["ci_high_95"] < 0
        else "negative"
        if result["ci_low_95"] > 0
        else "neutral"
    )
    result["material"] = abs(float(result["difference"])) >= 0.010
    return result


def repeat_comparison(predictions: pd.DataFrame) -> dict[str, object]:
    rows = []
    for method in ("A_structure_only", "F_full_solvai"):
        subset = predictions.loc[predictions.method.eq(method)].copy()
        subset["absolute_error"] = np.abs(subset.y_true - subset.y_pred)
        average = subset.groupby("molecule_id", as_index=False).agg(
            y_true=("y_true", "first"), absolute_error=("absolute_error", "mean")
        )
        average["method"] = method
        rows.append(average)
    combined = pd.concat(rows)
    candidate = combined.loc[combined.method.eq("F_full_solvai")].sort_values("molecule_id")
    baseline = combined.loc[combined.method.eq("A_structure_only")].sort_values("molecule_id")
    result: dict[str, object] = {
        "analysis": "five_repeat_molecule_averaged",
        "candidate": "F_full_solvai",
        "reference": "A_structure_only",
        "candidate_mae": float(candidate.absolute_error.mean()),
        "reference_mae": float(baseline.absolute_error.mean()),
    }
    result.update(
        paired_bootstrap(candidate.absolute_error.to_numpy() - baseline.absolute_error.to_numpy())
    )
    result["outcome"] = (
        "positive"
        if result["ci_high_95"] < 0
        else "negative"
        if result["ci_low_95"] > 0
        else "neutral"
    )
    result["material"] = abs(float(result["difference"])) >= 0.010
    return result


def shuffle_comparisons(aligned: pd.DataFrame, shuffled: pd.DataFrame) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    partitions = [("primary", -1)] + [("repeat", repeat) for repeat in range(5)]
    for partition, repeat in partitions:
        aligned_partition = aligned.loc[
            aligned.partition.eq(f"standardized_exclusion_{partition}")
            & aligned.repeat.eq(repeat)
            & aligned.method.eq("F_full_solvai")
        ].sort_values("molecule_id")
        shuffled_partition = shuffled.loc[
            shuffled.partition.eq(partition) & shuffled.repeat.eq(repeat)
        ].copy()
        shuffled_mean = (
            shuffled_partition.groupby("molecule_id", as_index=False)
            .agg(y_true=("y_true", "first"), y_pred=("y_pred", "mean"))
            .sort_values("molecule_id")
        )
        if len(aligned_partition) != 85 or len(shuffled_mean) != 85:
            raise AssertionError(f"Incomplete shuffled partition {partition}/{repeat}")
        aligned_error = np.abs(
            aligned_partition.y_true.to_numpy() - aligned_partition.y_pred.to_numpy()
        )
        shuffled_error = np.abs(shuffled_mean.y_true.to_numpy() - shuffled_mean.y_pred.to_numpy())
        result: dict[str, object] = {
            "analysis": f"aligned_vs_mean_shuffle_{partition}_{repeat}",
            "candidate": "F_full_solvai",
            "reference": "mean_shuffled_full_priors",
            "candidate_mae": float(aligned_error.mean()),
            "reference_mae": float(shuffled_error.mean()),
        }
        result.update(paired_bootstrap(aligned_error - shuffled_error))
        result["outcome"] = (
            "positive"
            if result["ci_high_95"] < 0
            else "negative"
            if result["ci_low_95"] > 0
            else "neutral"
        )
        result["material"] = abs(float(result["difference"])) >= 0.010
        results.append(result)
    return results


def format_ci(row: pd.Series | dict[str, object]) -> str:
    difference_key = "difference" if "difference" in row else "mae_difference_vs_structure"
    return (
        f"{float(row[difference_key]):.3f} "
        f"[{float(row['ci_low_95']):.3f}, {float(row['ci_high_95']):.3f}]"
    )


def main() -> None:
    endpoint_predictions = pd.read_parquet(
        OUT / "standardized_exclusion_endpoint_predictions.parquet"
    )
    shuffle_predictions = pd.read_parquet(
        OUT / "standardized_exclusion_endpoint_shuffle_predictions.parquet"
    )
    block_metrics = pd.read_csv(OUT / "standardized_exclusion_endpoint_metrics.csv")
    separation_metrics = pd.read_csv(OUT / "standardized_exclusion_global_separation_metrics.csv")
    separation_comparisons = pd.read_csv(
        OUT / "standardized_exclusion_global_separation_paired_comparisons.csv"
    )

    primary = endpoint_predictions.loc[
        endpoint_predictions.partition.eq("standardized_exclusion_primary")
    ]
    comparisons: list[dict[str, object]] = []
    for method in sorted(primary.method.unique()):
        if method == "A_structure_only":
            continue
        comparisons.append(
            comparison_row(
                f"primary_{method}",
                primary,
                primary,
                method,
                "A_structure_only",
            )
        )

    repeated = endpoint_predictions.loc[
        endpoint_predictions.partition.eq("standardized_exclusion_repeat")
    ]
    comparisons.append(repeat_comparison(repeated))
    comparisons.extend(shuffle_comparisons(endpoint_predictions, shuffle_predictions))

    zero_arrow = endpoint_predictions.loc[
        endpoint_predictions.partition.eq("standardized_exclusion_zero_arrow")
    ]
    comparisons.append(
        comparison_row(
            "zero_arrow_transfer",
            zero_arrow,
            zero_arrow,
            "F_full_solvai",
            "A_structure_only",
        )
    )
    comparison_table = pd.DataFrame(comparisons)
    comparison_table.to_csv(OUT / "confirmatory_paired_comparisons.csv", index=False)

    repeat_metrics = block_metrics.loc[
        block_metrics.partition.eq("standardized_exclusion_repeat")
    ].copy()
    repeat_summary = (
        repeat_metrics.groupby("method", as_index=False)
        .agg(
            repeats=("mae", "size"),
            mean_mae=("mae", "mean"),
            sd_mae=("mae", "std"),
            minimum_mae=("mae", "min"),
            maximum_mae=("mae", "max"),
        )
        .sort_values("method")
    )
    repeat_summary.to_csv(OUT / "confirmatory_repeat_summary.csv", index=False)

    shuffle_metrics = pd.read_csv(OUT / "standardized_exclusion_endpoint_shuffle_metrics.csv")
    shuffle_summary = (
        shuffle_metrics.groupby(["partition", "repeat"], as_index=False)
        .agg(mean_mae=("mae", "mean"), sd_mae=("mae", "std"))
        .sort_values(["partition", "repeat"])
    )
    shuffle_summary.to_csv(OUT / "confirmatory_shuffle_summary.csv", index=False)

    old_paper = json.loads((ROOT / "results" / "paper_metrics.json").read_text())
    pimd_mae = float(old_paper["methods"]["arrow_pimd8"]["mae_kcal_mol"])
    primary_metrics = block_metrics.loc[
        block_metrics.partition.eq("standardized_exclusion_primary")
    ].set_index("method")
    full_primary = float(primary_metrics.loc["F_full_solvai", "mae"])
    structure_primary = float(primary_metrics.loc["A_structure_only", "mae"])
    full_repeats = repeat_summary.set_index("method").loc["F_full_solvai"]
    structure_repeats = repeat_summary.set_index("method").loc["A_structure_only"]

    headline = {
        "schema_version": 1,
        "protocol": "release/CONFIRMATORY_FREEZE.md",
        "benchmark_molecules": 85,
        "external_endpoint_labels": 1280,
        "standardized_teacher_exclusions": {
            "combisolv_qm": 2,
            "molsolv_smd": 32,
            "confsolv": 22,
        },
        "primary": {
            "matched_structure_only_mae": structure_primary,
            "full_solvai_mae": full_primary,
            "mae_difference": full_primary - structure_primary,
            "arrow_pimd8_mae": pimd_mae,
        },
        "repeats": {
            "structure_only_values": repeat_metrics.loc[
                repeat_metrics.method.eq("A_structure_only"), "mae"
            ].tolist(),
            "structure_only_mean": float(structure_repeats.mean_mae),
            "structure_only_sd": float(structure_repeats.sd_mae),
            "full_solvai_values": repeat_metrics.loc[
                repeat_metrics.method.eq("F_full_solvai"), "mae"
            ].tolist(),
            "full_solvai_mean": float(full_repeats.mean_mae),
            "full_solvai_sd": float(full_repeats.sd_mae),
        },
        "zero_arrow_transfer": {
            row.method: float(row.mae)
            for row in block_metrics.loc[
                block_metrics.partition.eq("standardized_exclusion_zero_arrow")
            ].itertuples()
        },
        "global_separation": separation_metrics.to_dict(orient="records"),
        "historical_context_not_confirmatory": {
            "reported_structure_only_label_value": float(
                old_paper["methods"]["previous_structure_only"]["mae_kcal_mol"]
            ),
            "reported_fixed_smd_confsolv_value": float(
                old_paper["methods"]["smd_confsolv_fixed"]["mae_kcal_mol"]
            ),
            "note": (
                "The historical 0.238606 comparator was a physics-mixture model, "
                "not the matched descriptor-only endpoint required by Phase 1."
            ),
        },
        "optional_analyses": {
            "learning_curve": "not run; core fixed analyses took priority",
            "teacher_fidelity_correlation": (
                "not run as a cross-source causal analysis because targets, units "
                "and test protocols are not commensurate"
            ),
        },
    }
    (OUT / "confirmatory_summary.json").write_text(json.dumps(headline, indent=2) + "\n")

    summary_rows = [
        {
            "analysis": "primary_fixed",
            "method": row.method,
            "mae": row.mae,
            "rmse": row.rmse,
            "n": row.n,
        }
        for row in primary_metrics.reset_index().itertuples()
    ]
    summary_rows.extend(
        {
            "analysis": row.regime,
            "method": row.method,
            "mae": row.mae,
            "rmse": row.rmse,
            "n": row.n,
        }
        for row in separation_metrics.itertuples()
    )
    pd.DataFrame(summary_rows).to_csv(OUT / "confirmatory_summary.csv", index=False)

    primary_full_comparison = comparison_table.loc[
        comparison_table.analysis.eq("primary_F_full_solvai")
    ].iloc[0]
    repeat_comparison_row = comparison_table.loc[
        comparison_table.analysis.eq("five_repeat_molecule_averaged")
    ].iloc[0]
    zero_comparison = comparison_table.loc[
        comparison_table.analysis.eq("zero_arrow_transfer")
    ].iloc[0]
    primary_shuffle = comparison_table.loc[
        comparison_table.analysis.eq("aligned_vs_mean_shuffle_primary_-1")
    ].iloc[0]

    block_lines = []
    for row in primary_metrics.reset_index().itertuples():
        block_lines.append(f"| {row.method} | {row.mae:.3f} | {row.rmse:.3f} |")
    separation_lines = []
    separated = separation_metrics.pivot(index="regime", columns="method", values="mae")
    paired_sep = separation_comparisons.set_index("regime")
    for regime, row in separated.iterrows():
        comparison = paired_sep.loc[regime]
        separation_lines.append(
            f"| {regime} | {row['A_structure_only']:.3f} | "
            f"{row['F_full_solvai']:.3f} | {format_ci(comparison)} |"
        )

    REPORT.write_text(
        f"""# SolvAI Phase 1 confirmatory analysis

## Executive result

The preregistered matched comparison supports the narrow scientific claim that
molecule-aligned, structure-predicted response priors add useful information to the
frozen endpoint pipeline. After removing standardized-equivalent benchmark records
from affected teacher sources and refitting those teachers with original split
membership preserved, the matched descriptor-only endpoint has an MAE of
**{structure_primary:.3f} kcal mol\u207b\u00b9**, whereas the full 15-prior SolvAI stack has an
MAE of **{full_primary:.3f} kcal mol\u207b\u00b9**. The paired difference is
{format_ci(primary_full_comparison)} kcal mol\u207b\u00b9 (candidate minus baseline;
95% molecule bootstrap interval).

This is a stronger causal control than the prior campaign comparison. It also
changes the publication record in two important ways: the historical 0.238606 value
was not a matched descriptor-only endpoint, and the historical 0.197047 model used
teachers filtered only by exact connectivity. The standardized-exclusion
confirmatory point estimate, **{full_primary:.3f}**, is the appropriate primary
fixed-partition value. It is on the ARROW/PIMD8 accuracy scale
({pimd_mae:.3f}) but is not evidence of superiority.

## Predeclared endpoint controls

| Fixed feature set | MAE | RMSE |
| --- | ---: | ---: |
{chr(10).join(block_lines)}

The blocks are scientific ablations, not an additive attribution. In isolation,
the empirical/residual block improves the matched endpoint; SMD alone is positive
under the original exact-connectivity teachers but becomes statistically neutral
after the more conservative standardized exclusion. ConfSolv alone is neutral. The
full block is positive, consistent with complementary information across response
coordinates rather than one universally dominant teacher.

Across the five preregistered partitions, the descriptor-only endpoint is
**{float(structure_repeats.mean_mae):.3f} \u00b1 {float(structure_repeats.sd_mae):.3f}** and
full SolvAI is **{float(full_repeats.mean_mae):.3f} \u00b1 {float(full_repeats.sd_mae):.3f}**
(mean \u00b1 sample s.d.). The molecule-averaged paired improvement is
{format_ci(repeat_comparison_row)} kcal mol\u207b\u00b9. Full SolvAI improves on the matched
endpoint in all five partitions.

## Shuffled-prior control

Meaningful molecule--response alignment is necessary. On the fixed partition,
the mean of five shuffled-prior controls is **{float(primary_shuffle.reference_mae):.3f}**
MAE versus **{float(primary_shuffle.candidate_mae):.3f}** for aligned priors; the
paired difference is {format_ci(primary_shuffle)} kcal mol\u207b\u00b9. Shuffled priors
remain near the matched structure-only endpoint rather than reproducing the aligned
gain. This rejects the explanation that arbitrary extra columns or model width are
sufficient.

## Global chemical separation

Chemical separation was applied to every endpoint-supervised molecule, including
the 1,280-label external pool.

| Regime | Structure-only MAE | Full SolvAI MAE | Paired difference [95% CI] |
| --- | ---: | ---: | ---: |
{chr(10).join(separation_lines)}

The prior advantage survives every predeclared separation regime, with paired
intervals below zero. Absolute error rises sharply for leave-family and scaffold
extrapolation, so these results support transfer of the response-prior advantage,
not broad high-accuracy extrapolation.

## Zero-ARROW-label transfer

With no ARROW experimental label used in endpoint fitting, the matched structure
model gives **{float(zero_comparison.reference_mae):.3f}** MAE and full SolvAI gives
**{float(zero_comparison.candidate_mae):.3f}**. The paired difference is
{format_ci(zero_comparison)} kcal mol\u207b\u00b9. This supports representation value without
fold-local adaptation, but it is not an independent external validation because
the method was developed in the context of the 85-solute set.

## Identity and chemical-distance audit

The endpoint pool contains no exact, fragment-parent, uncharged-parent or canonical-
tautomer match to the 85 reference connectivities. The expanded audit found
standardized equivalents in three teacher sources: 2 CombiSolv-QM rows, 32 MolSolv
records and 22 ConfSolv rows. Those records were excluded, the affected teachers
were refitted while preserving every remaining molecule's original source split,
and all affected endpoint analyses were repeated. Morgan-similarity 1.0 collisions
between non-identical structures are reported as fingerprint collisions, not hidden
identity matches. Full details and every pair appear under `audits/confirmatory/`.

## Claims that survived

- The 15 aligned response priors materially improve an otherwise matched endpoint.
- Arbitrarily shuffled priors do not reproduce the gain.
- The advantage survives family, scaffold, cluster and nearest-neighbour separation
  applied globally to all endpoint labels.
- The advantage remains when no ARROW label is used for endpoint fitting.
- Structure-only SolvAI reaches the ARROW/PIMD8 accuracy scale on the fixed reference
  partition, with no simulation at inference.

## Claims that weakened or died

- The historical 0.238606 number cannot be described as the matched no-prior
  baseline; the correct matched value is {structure_primary:.3f}.
- The historical 0.197047 estimate is not the conservative publication headline
  after the expanded identity audit; the corrected fixed result is {full_primary:.3f}.
- Robust sub-0.20 performance is not supported. Corrected repeated performance
  centres at {float(full_repeats.mean_mae):.3f} kcal mol\u207b\u00b9.
- Chemical extrapolation is not solved: global family and scaffold MAEs are
  {float(separated.loc["global_family", "F_full_solvai"]):.3f} and
  {float(separated.loc["global_scaffold", "F_full_solvai"]):.3f}, respectively.
- PIMD supervision was not retained; the result concerns distillation of diverse
  solvation-response calculations, with PIMD8 serving only as the accuracy
  comparator.

## Nature Communications thesis

The confirmatory evidence supports a focused thesis: external physical calculations
can define molecule-aligned response coordinates that are learned from structure and
provide a reproducible endpoint advantage beyond the same experimental labels,
representation, model class, folds and seeds. It does not establish a universal
principle across properties or chemistry. The manuscript should therefore lead with
the controlled response-prior result and use the PIMD8 comparison as a scientifically
important scale reference, not as the causal proof.

## Execution record

All definitions are frozen in `release/CONFIRMATORY_FREEZE.md`. Valid analyses used
`uv run` with Python 3.11.15, RDKit 2026.03.5 and scikit-learn 1.7.2. Primary commands:

```bash
uv run python scripts/run_confirmatory_endpoint.py --mode primary
uv run python scripts/run_confirmatory_endpoint.py --mode repeats
uv run python scripts/run_confirmatory_endpoint.py --mode shuffle
uv run python scripts/audit_confirmatory_chemistry.py
uv run python scripts/train_standardized_exclusion_teachers.py
uv run python scripts/verify_confirmatory_teacher_refits.py
uv run python scripts/run_standardized_exclusion_endpoints.py
uv run python scripts/run_confirmatory_endpoint.py --mode shuffle --standardized-exclusion
uv run python scripts/run_confirmatory_separation.py --standardized-exclusion
uv run python scripts/summarize_confirmatory.py
```

The optional experimental-label learning curve was not run because it required 750
additional endpoint fits and was not needed to resolve any primary interpretation
rule. A cross-source teacher-fidelity regression was not run because the available
teacher targets differ in units, scales, model classes and test protocols; treating
them as commensurate would be more misleading than informative.
"""
    )


if __name__ == "__main__":
    main()
