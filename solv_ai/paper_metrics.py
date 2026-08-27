"""Canonical computation of every numerical result used by the SolvAI paper.

The manuscript, tables and figures consume ``results/paper_metrics.json`` and
``results/paper_metrics.csv`` generated here.  This module reads only frozen,
molecule-level outputs and asserts their identities, row counts and stored-error
consistency before writing publication values.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP_REPLICATES = 100_000
BOOTSTRAP_SEED = 20260827

METHOD = {
    "previous": "Nested structure-only physics mixture",
    "response": "Fixed narrow response without SMD",
    "smd": "Fixed narrow response + SMD water",
    "final": "Fixed narrow response + SMD + ConfSolv response",
    "nested": "Nested narrow SMD teacher selection",
    "hard": "SMD response direct + ConfSolv response",
    "pimd_residual": "PIMD8 + ML residual (nested)",
}


def _prediction_block(
    frame: pd.DataFrame, method: str, regime: str, *, repeats: bool = False
) -> pd.DataFrame:
    block = frame.loc[frame.method.eq(method) & frame.regime.eq(regime)].copy()
    expected = 425 if repeats else 85
    if len(block) != expected:
        raise AssertionError(f"{method}/{regime}: expected {expected} rows, found {len(block)}")
    key = ["repeat", "molecule_id"] if repeats else ["molecule_id"]
    if block.duplicated(key).any():
        raise AssertionError(f"{method}/{regime}: duplicate prediction identities")
    recomputed = (block.y_pred.astype(float) - block.y_true.astype(float)).abs()
    if not np.allclose(recomputed, block.absolute_error.astype(float), atol=1e-12):
        raise AssertionError(f"{method}/{regime}: stale absolute-error column")
    return block


def _summary(block: pd.DataFrame) -> dict[str, float | int]:
    error = block.y_pred.to_numpy(float) - block.y_true.to_numpy(float)
    return {
        "n": len(block),
        "mae_kcal_mol": float(np.mean(np.abs(error))),
        "rmse_kcal_mol": float(np.sqrt(np.mean(np.square(error)))),
    }


def _bootstrap(errors: np.ndarray, *, seed: int = BOOTSTRAP_SEED) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    distribution = np.empty(BOOTSTRAP_REPLICATES, dtype=float)
    for start in range(0, BOOTSTRAP_REPLICATES, 10_000):
        stop = min(start + 10_000, BOOTSTRAP_REPLICATES)
        index = rng.integers(0, len(errors), size=(stop - start, len(errors)))
        distribution[start:stop] = errors[index].mean(axis=1)
    low, high = np.quantile(distribution, [0.025, 0.975])
    return {
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": seed,
        "ci95_kcal_mol": [float(low), float(high)],
        "probability_mae_below_0_20": float(np.mean(distribution < 0.20)),
    }


def _paired_bootstrap(
    candidate_error: np.ndarray, reference_error: np.ndarray, *, seed: int
) -> dict[str, Any]:
    if candidate_error.shape != reference_error.shape:
        raise AssertionError("Paired errors must have identical shapes")
    delta = candidate_error - reference_error
    rng = np.random.default_rng(seed)
    distribution = np.empty(BOOTSTRAP_REPLICATES, dtype=float)
    for start in range(0, BOOTSTRAP_REPLICATES, 10_000):
        stop = min(start + 10_000, BOOTSTRAP_REPLICATES)
        index = rng.integers(0, len(delta), size=(stop - start, len(delta)))
        distribution[start:stop] = delta[index].mean(axis=1)
    low, high = np.quantile(distribution, [0.025, 0.975])
    return {
        "mean_mae_change_kcal_mol": float(delta.mean()),
        "ci95_kcal_mol": [float(low), float(high)],
        "probability_candidate_better": float(np.mean(distribution < 0.0)),
    }


def compute_paper_metrics(root: Path = ROOT) -> tuple[dict[str, Any], pd.DataFrame]:
    benchmark = pd.read_parquet(root / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")].sort_values("molecule_id")
    if len(benchmark) != 85 or benchmark.molecule_id.duplicated().any():
        raise AssertionError("The frozen ARROW hydration reference must contain 85 molecules")

    headline = pd.read_parquet(root / "results/predictions/headline_oof.parquet")
    hard = pd.read_parquet(root / "results/predictions/hard_holdout_oof.parquet")
    repeated = pd.read_parquet(root / "results/robustness/repeated_oof.parquet")
    campaign = pd.read_parquet(root / "results/predictions/campaign_oof.parquet")

    blocks = {
        key: _prediction_block(headline, METHOD[key], "random_oof")
        for key in ("response", "smd", "final", "nested")
    }
    previous = _prediction_block(campaign, METHOD["previous"], "random_oof")
    pimd_residual = _prediction_block(campaign, METHOD["pimd_residual"], "random_oof")
    family = _prediction_block(hard, METHOD["hard"], "family_holdout")
    scaffold = _prediction_block(hard, METHOD["hard"], "scaffold_holdout")

    ids = benchmark.molecule_id.astype(str).tolist()
    for name, block in {**blocks, "previous": previous, "pimd_residual": pimd_residual}.items():
        if sorted(block.molecule_id.astype(str)) != ids:
            raise AssertionError(f"{name}: OOF identities differ from frozen benchmark")

    truth = benchmark.delta_g_exp.to_numpy(float)
    pimd_error = benchmark.delta_g_pimd8.to_numpy(float) - truth
    classical_error = benchmark.delta_g_classical_arrow.to_numpy(float) - truth
    method_metrics = {
        "classical_arrow": {
            "n": 85,
            "mae_kcal_mol": float(np.mean(np.abs(classical_error))),
            "rmse_kcal_mol": float(np.sqrt(np.mean(np.square(classical_error)))),
        },
        "arrow_pimd8": {
            "n": 85,
            "mae_kcal_mol": float(np.mean(np.abs(pimd_error))),
            "rmse_kcal_mol": float(np.sqrt(np.mean(np.square(pimd_error)))),
        },
        "pimd8_plus_nested_residual": _summary(pimd_residual),
        "previous_structure_only": _summary(previous),
        "narrow_response": _summary(blocks["response"]),
        "smd_water": _summary(blocks["smd"]),
        "smd_confsolv_fixed": _summary(blocks["final"]),
        "nested_selection": _summary(blocks["nested"]),
        "family_holdout": _summary(family),
        "scaffold_holdout": _summary(scaffold),
    }

    repeat_metrics: dict[str, Any] = {}
    for key, method in {
        "fixed": "narrow response + SMD + ConfSolv response",
        "nested": "Nested selection",
    }.items():
        block = repeated.loc[repeated.method.eq(method)].copy()
        if len(block) != 425 or block.duplicated(["repeat", "molecule_id"]).any():
            raise AssertionError(f"{method}: expected five complete 85-molecule repeats")
        per_repeat = block.groupby("repeat", sort=True).absolute_error.mean()
        if per_repeat.index.tolist() != [0, 1, 2, 3, 4]:
            raise AssertionError(f"{method}: repeat IDs are incomplete")
        repeat_metrics[key] = {
            "values_kcal_mol": [float(value) for value in per_repeat],
            "mean_kcal_mol": float(per_repeat.mean()),
            "sd_kcal_mol": float(per_repeat.std(ddof=1)),
            "repeats_below_0_20": int(per_repeat.lt(0.20).sum()),
        }

    fixed = blocks["final"].sort_values("molecule_id")
    nested = blocks["nested"].sort_values("molecule_id")
    prior = previous.sort_values("molecule_id")
    pimd_abs = np.abs(pimd_error)
    fixed_abs = fixed.absolute_error.to_numpy(float)
    nested_abs = nested.absolute_error.to_numpy(float)
    prior_abs = prior.absolute_error.to_numpy(float)
    uncertainty = {
        "fixed": _bootstrap(fixed_abs, seed=BOOTSTRAP_SEED),
        "nested": _bootstrap(nested_abs, seed=BOOTSTRAP_SEED + 1),
        "fixed_vs_previous": _paired_bootstrap(fixed_abs, prior_abs, seed=BOOTSTRAP_SEED + 2),
        "fixed_vs_pimd8": _paired_bootstrap(fixed_abs, pimd_abs, seed=BOOTSTRAP_SEED + 3),
    }

    fixed_family = (
        fixed.groupby("functional_group_family", sort=True)
        .absolute_error.agg(["count", "mean"])
        .reset_index()
    )
    chemistry_family = [
        {
            "family": str(row.functional_group_family),
            "n": int(row["count"]),
            "mae_kcal_mol": float(row["mean"]),
        }
        for _, row in fixed_family.iterrows()
    ]

    multilambda = pd.read_csv(root / "results/ablations/multilambda_metrics.csv")
    response_metrics = pd.read_csv(root / "results/ablations/multilambda_response_metrics.csv")
    lambda_teachers = pd.read_parquet(root / "results/ablations/pimd2_multilambda_teacher.parquet")
    lambda_counts = {
        code: int(lambda_teachers[f"success__lambda{code}"].fillna(False).sum())
        for code in ("01", "05", "09")
    }
    lambda_counts["complete_three_point"] = int(lambda_teachers.complete_three_point_curve.sum())
    lambda_results = {
        str(row.method): float(row.mae) for row in multilambda.itertuples(index=False)
    }
    response_records = response_metrics.replace({np.nan: None}).to_dict(orient="records")

    screens = {}
    for key, filename in {
        "matched_one_seed_base": "smd_plus_confsolv_response_screen.csv",
        "openfe_diagnostics": "openfe_screen.csv",
        "mlff_hierarchy": "mlff_screen.csv",
        "des370k_water_response": "des370k_screen.csv",
    }.items():
        screen = pd.read_csv(root / "results/ablations" / filename)
        screens[key] = float(screen.mae.iloc[0])
    selected_alternatives = {
        "qmpff_typed_structure": (
            "qmpff_static_ablation_oof.csv",
            "expanded repeated/curated; legacy + Abraham + QMPFF typed structure",
        ),
        "gnnis_force_response": (
            "smd_plus_gnnis_force_structure_screen.csv",
            "SMD response direct + GNNIS force",
        ),
        "lambda_aware_implicit": (
            "smd_plus_lambda_potential_structure_screen.csv",
            "SMD response direct + lambda potential",
        ),
        "phase_space_dynamic": (
            "new_phase_hierarchy_screen_oof.csv",
            "expanded repeated/curated; legacy + all response + phase-space dynamic",
        ),
        "gnequip_hierarchy": (
            "new_phase_hierarchy_screen_oof.csv",
            "expanded repeated/curated; legacy + all response + G-NequIP hierarchy",
        ),
        "confsolv_graph_embedding": (
            "smd_confsolv_graph_embedding_screen.csv",
            "SMD response direct + ConfSolv response + ConfSolv graph embedding",
        ),
        "confsolv_ffn_embedding": (
            "smd_confsolv_ffn_embedding_screen.csv",
            "SMD response direct + ConfSolv response + ConfSolv ffn embedding",
        ),
        "openfe_nqe_dpnn": (
            "foundation_openfe_nqe_hierarchy_exp_pimd_oof_predictions.csv",
            "D-MPNN openfe_nqe_hierarchy exp_pimd",
        ),
        "molpile_dpnn": (
            "foundation_molpile10m_exp_oof_predictions.csv",
            "D-MPNN molpile10m exp",
        ),
    }
    for key, (filename, method) in selected_alternatives.items():
        table = pd.read_csv(root / "results/ablations" / filename)
        row = table.loc[table.method.eq(method) & table.regime.eq("random_oof")]
        if len(row) != 1:
            raise AssertionError(f"Missing frozen alternative result: {key}")
        screens[key] = float(row.mae.iloc[0] if "mae" in row else row["mean"].iloc[0])

    molsolv_metadata = json.loads(
        (root / "data/manifests/molsolv_smd_teacher_metadata.json").read_text()
    )
    confsolv_metadata = json.loads(
        (root / "data/manifests/confsolv_water_teacher_metadata.json").read_text()
    )
    molsolv = pd.read_parquet(root / "data/processed/molsolv_smd_water_nonbenchmark.parquet")
    confsolv = pd.read_parquet(root / "data/processed/confsolv_water_nonbenchmark.parquet")
    endpoint_metadata = json.loads(
        (root / "data/manifests/endpoint_label_manifest.json").read_text()
    )
    data_counts = {
        "molsolv_source_conformers": int(molsolv_metadata["records_total"]),
        "molsolv_training_structures": len(molsolv),
        "molsolv_unique_connectivities": int(molsolv.connectivity_key.nunique()),
        "molsolv_benchmark_structure_matches_removed": int(
            molsolv_metadata["benchmark_exact_structure_matches_removed"]
        ),
        "confsolv_source_water_conformers": int(confsolv_metadata["raw_water_rows"]),
        "confsolv_training_connectivities": len(confsolv),
        "confsolv_benchmark_connectivities_removed": int(
            confsolv_metadata["benchmark_connectivity_overlaps_removed"]
        ),
        "public_hydration_connectivities": int(endpoint_metadata["expanded_table"]["rows"]),
        "endpoint_public_labels": int(endpoint_metadata["endpoint_selection"]["rows"]),
        "arrow85_connectivities_also_in_freesolv": int(benchmark.freesolv_match_id.notna().sum()),
    }

    artifact = joblib.load(root / "models/final/head.joblib")
    artifact_metrics = {
        "input": "SMILES only",
        "simulation_at_inference": False,
        "descriptor_features": len(artifact["descriptor_columns"]),
        "physics_response_features": len(artifact["teacher_columns"]),
        "total_features": int(
            len(artifact["descriptor_columns"]) + len(artifact["teacher_columns"])
        ),
        "ensemble_members": len(artifact["models"]),
        "seeds": [int(seed) for seed in artifact["seeds"]],
    }

    expected = {
        "arrow_pimd8": 0.20483647058823526,
        "previous_structure_only": 0.23860611898039194,
        "narrow_response": 0.2136151731870528,
        "smd_water": 0.2016070042016807,
        "smd_confsolv_fixed": 0.197047474094823,
        "nested_selection": 0.1993060393064633,
        "family_holdout": 0.2395665848765437,
        "scaffold_holdout": 0.2412820403309475,
    }
    for key, value in expected.items():
        observed = method_metrics[key]["mae_kcal_mol"]
        if not np.isclose(observed, value, atol=1e-12):
            raise AssertionError(f"Frozen metric drift for {key}: {observed} != {value}")
    if data_counts["molsolv_training_structures"] != 350_391:
        raise AssertionError("MolSolv filtered count drift")
    if data_counts["confsolv_training_connectivities"] != 39_878:
        raise AssertionError("ConfSolv filtered count drift")
    if data_counts["public_hydration_connectivities"] != 5_075:
        raise AssertionError("Public hydration catalog count drift")
    if data_counts["endpoint_public_labels"] != 1_280:
        raise AssertionError("Endpoint label-selection count drift")
    if data_counts["arrow85_connectivities_also_in_freesolv"] != 80:
        raise AssertionError("ARROW/FreeSolv identity-overlap count drift")
    if artifact_metrics["total_features"] != 2_280:
        raise AssertionError("Packaged feature schema drift")

    metrics: dict[str, Any] = {
        "schema_version": 1,
        "benchmark": {
            "name": "ARROW 85-solute neutral-hydration reference set",
            "molecules": 85,
            "temperature_k": 298.0,
            "freesolv_connectivity_overlap": 80,
        },
        "methods": method_metrics,
        "repeated_splits": repeat_metrics,
        "bootstrap": uncertainty,
        "chemistry_family": chemistry_family,
        "alternative_supervision": screens,
        "multilambda": {
            "teacher_counts": lambda_counts,
            "method_mae_kcal_mol": lambda_results,
            "response_head_metrics": response_records,
        },
        "data_counts": data_counts,
        "artifact": artifact_metrics,
        "audit": {
            "assertions_passed": True,
            "metric_tolerance": 1e-12,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "bootstrap_seed": BOOTSTRAP_SEED,
        },
    }

    tidy_rows: list[dict[str, Any]] = []
    for key, value in method_metrics.items():
        tidy_rows.extend(
            [
                {
                    "section": "method",
                    "metric": f"{key}.mae",
                    "value": value["mae_kcal_mol"],
                    "unit": "kcal/mol",
                    "n": value["n"],
                },
                {
                    "section": "method",
                    "metric": f"{key}.rmse",
                    "value": value["rmse_kcal_mol"],
                    "unit": "kcal/mol",
                    "n": value["n"],
                },
            ]
        )
    for key, value in repeat_metrics.items():
        tidy_rows.extend(
            [
                {
                    "section": "repeat",
                    "metric": f"{key}.mean_mae",
                    "value": value["mean_kcal_mol"],
                    "unit": "kcal/mol",
                    "n": 5,
                },
                {
                    "section": "repeat",
                    "metric": f"{key}.sd_mae",
                    "value": value["sd_kcal_mol"],
                    "unit": "kcal/mol",
                    "n": 5,
                },
            ]
        )
        for index, score in enumerate(value["values_kcal_mol"]):
            tidy_rows.append(
                {
                    "section": "repeat",
                    "metric": f"{key}.repeat_{index}",
                    "value": score,
                    "unit": "kcal/mol",
                    "n": 85,
                }
            )
    for key, value in data_counts.items():
        tidy_rows.append(
            {"section": "data", "metric": key, "value": value, "unit": "count", "n": np.nan}
        )
    for key, value in artifact_metrics.items():
        if isinstance(value, (int, float, bool)):
            tidy_rows.append(
                {
                    "section": "artifact",
                    "metric": key,
                    "value": value,
                    "unit": "count" if isinstance(value, int) else "",
                    "n": np.nan,
                }
            )
    return metrics, pd.DataFrame(tidy_rows)


def freeze_markdown(metrics: dict[str, Any]) -> str:
    methods = metrics["methods"]
    repeats = metrics["repeated_splits"]
    data = metrics["data_counts"]
    return f"""# SolvAI paper freeze

This file freezes the quantitative state used by the manuscript. It is generated
only by `solv_ai.paper_metrics`; hand-edited numerical values are not authoritative.

## Scientific decision

SolvAI is a SMILES-only inference stack. The strict fixed five-fold OOF result is
**{methods["smd_confsolv_fixed"]["mae_kcal_mol"]:.5f} kcal/mol**, and nested
feature-block selection is **{methods["nested_selection"]["mae_kcal_mol"]:.5f}**.
The fixed model averages **{repeats["fixed"]["mean_kcal_mol"]:.5f} ±
{repeats["fixed"]["sd_kcal_mol"]:.5f}** across five independent split repeats.
Accordingly, the single-partition threshold crossing is valid but robust sub-0.20
generalization is not claimed.

## Frozen values

| Quantity | Value |
|---|---:|
| Molecules | {metrics["benchmark"]["molecules"]} |
| Classical ARROW MAE | {methods["classical_arrow"]["mae_kcal_mol"]:.5f} kcal/mol |
| ARROW/PIMD8 MAE | {methods["arrow_pimd8"]["mae_kcal_mol"]:.5f} kcal/mol |
| Previous structure-only MAE | {methods["previous_structure_only"]["mae_kcal_mol"]:.5f} kcal/mol |
| Narrow-response MAE | {methods["narrow_response"]["mae_kcal_mol"]:.5f} kcal/mol |
| + MolSolv SMD(water) MAE | {methods["smd_water"]["mae_kcal_mol"]:.5f} kcal/mol |
| + ConfSolv response MAE | {methods["smd_confsolv_fixed"]["mae_kcal_mol"]:.5f} kcal/mol |
| Nested-selection MAE | {methods["nested_selection"]["mae_kcal_mol"]:.5f} kcal/mol |
| Five-repeat fixed mean ± SD | {repeats["fixed"]["mean_kcal_mol"]:.5f} ± {repeats["fixed"]["sd_kcal_mol"]:.5f} kcal/mol |
| Five-repeat nested mean ± SD | {repeats["nested"]["mean_kcal_mol"]:.5f} ± {repeats["nested"]["sd_kcal_mol"]:.5f} kcal/mol |
| Family-held-out MAE | {methods["family_holdout"]["mae_kcal_mol"]:.5f} kcal/mol |
| Scaffold-held-out MAE | {methods["scaffold_holdout"]["mae_kcal_mol"]:.5f} kcal/mol |
| MolSolv source / retained | {data["molsolv_source_conformers"]:,} / {data["molsolv_training_structures"]:,} |
| ConfSolv source / retained | {data["confsolv_source_water_conformers"]:,} / {data["confsolv_training_connectivities"]:,} |
| ARROW-85 connectivities also in FreeSolv | {data["arrow85_connectivities_also_in_freesolv"]}/85 |
| Simulation at inference | No |

## Integrity

- Every OOF method has exactly 85 unique molecule predictions.
- Every repeated method has five complete 85-molecule partitions.
- Stored absolute errors equal values recomputed from predictions and targets.
- External count and model-schema assertions pass before files are written.
- Bootstrap interval: molecule-level resampling, {BOOTSTRAP_REPLICATES:,} replicates,
  seed {BOOTSTRAP_SEED}.
"""


def write_paper_metrics(root: Path = ROOT) -> dict[str, Any]:
    metrics, tidy = compute_paper_metrics(root)
    results = root / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "paper_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    tidy.to_csv(results / "paper_metrics.csv", index=False)
    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "PAPER_FREEZE.md").write_text(freeze_markdown(metrics))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    metrics = write_paper_metrics(args.root.resolve())
    print(json.dumps(metrics["methods"], indent=2))


if __name__ == "__main__":
    main()
