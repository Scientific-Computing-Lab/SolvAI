"""Canonical computation of every numerical result used by the SolvAI paper.

The publication source consumes only ``results/paper_metrics.json`` and
``results/paper_metrics.csv`` generated here. Confirmatory metrics are recomputed
from molecule-level predictions produced under ``release/CONFIRMATORY_FREEZE.md``.
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
BOOTSTRAP_SEED = 20260828


def _summary(block: pd.DataFrame) -> dict[str, float | int]:
    error = block.y_pred.to_numpy(float) - block.y_true.to_numpy(float)
    return {
        "n": len(block),
        "mae_kcal_mol": float(np.mean(np.abs(error))),
        "rmse_kcal_mol": float(np.sqrt(np.mean(np.square(error)))),
        "median_absolute_error_kcal_mol": float(np.median(np.abs(error))),
    }


def _bootstrap(values: np.ndarray, *, seed: int) -> dict[str, Any]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    samples = np.empty(BOOTSTRAP_REPLICATES, dtype=float)
    for start in range(0, BOOTSTRAP_REPLICATES, 2_000):
        stop = min(start + 2_000, BOOTSTRAP_REPLICATES)
        index = rng.integers(0, len(values), size=(stop - start, len(values)))
        samples[start:stop] = values[index].mean(axis=1)
    low, high = np.quantile(samples, [0.025, 0.975])
    return {
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": seed,
        "mean": float(values.mean()),
        "ci95": [float(low), float(high)],
        "probability_below_zero": float(np.mean(samples < 0.0)),
    }


def _method_block(frame: pd.DataFrame, partition: str, method: str) -> pd.DataFrame:
    block = frame.loc[frame.partition.eq(partition) & frame.method.eq(method)].copy()
    if len(block) != 85 or block.molecule_id.duplicated().any():
        raise AssertionError(f"{partition}/{method}: expected 85 unique predictions")
    recomputed = np.abs(block.y_pred.astype(float) - block.y_true.astype(float))
    if not np.allclose(recomputed, block.absolute_error.astype(float), atol=1e-12):
        raise AssertionError(f"{partition}/{method}: stale error column")
    return block.sort_values("molecule_id").reset_index(drop=True)


def _legacy_method(
    frame: pd.DataFrame, method: str, regime: str = "random_oof"
) -> dict[str, float | int]:
    block = frame.loc[frame.method.eq(method) & frame.regime.eq(regime)].copy()
    if len(block) != 85:
        raise AssertionError(f"Legacy method {method}/{regime} is incomplete")
    return _summary(block)


def compute_paper_metrics(root: Path = ROOT) -> tuple[dict[str, Any], pd.DataFrame]:
    benchmark = pd.read_parquet(root / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")].sort_values("molecule_id")
    if len(benchmark) != 85 or benchmark.molecule_id.duplicated().any():
        raise AssertionError("The ARROW reference must contain 85 unique water rows")

    endpoint = pd.read_parquet(
        root / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    shuffle = pd.read_parquet(
        root / "results/confirmatory/standardized_exclusion_endpoint_shuffle_predictions.parquet"
    )
    separation = pd.read_parquet(
        root / "results/confirmatory/standardized_exclusion_global_separation_predictions.parquet"
    )
    paired = pd.read_csv(root / "results/confirmatory/confirmatory_paired_comparisons.csv")
    chemistry_audit = json.loads(
        (root / "audits/confirmatory/chemical_distance_audit.json").read_text()
    )

    primary_partition = "standardized_exclusion_primary"
    zero_partition = "standardized_exclusion_zero_arrow"
    method_names = {
        "matched_structure_only": "A_structure_only",
        "empirical_residual_block": "B_empirical_residual",
        "computation_core_block": "C_computation_core",
        "smd_water_block": "D_smd_water",
        "confsolv_block": "E_confsolv",
        "full_solvai": "F_full_solvai",
        "narrow_response": "G_narrow_reference",
        "narrow_plus_smd": "H_narrow_smd_reference",
    }
    blocks = {
        key: _method_block(endpoint, primary_partition, method)
        for key, method in method_names.items()
    }

    truth = benchmark.delta_g_exp.to_numpy(float)
    classical_error = benchmark.delta_g_classical_arrow.to_numpy(float) - truth
    pimd_error = benchmark.delta_g_pimd8.to_numpy(float) - truth
    methods: dict[str, dict[str, float | int]] = {
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
    }
    methods.update({key: _summary(block) for key, block in blocks.items()})
    methods["zero_arrow_structure_only"] = _summary(
        _method_block(endpoint, zero_partition, "A_structure_only")
    )
    methods["zero_arrow_full_solvai"] = _summary(
        _method_block(endpoint, zero_partition, "F_full_solvai")
    )

    separation_metrics: dict[str, dict[str, dict[str, float | int]]] = {}
    for regime in sorted(separation.regime.unique()):
        separation_metrics[regime] = {}
        for method in ("A_structure_only", "F_full_solvai"):
            block = separation.loc[
                separation.regime.eq(regime) & separation.method.eq(method)
            ].copy()
            if len(block) != 85:
                raise AssertionError(f"{regime}/{method}: incomplete separation result")
            separation_metrics[regime][method] = _summary(block)

    repeated = endpoint.loc[endpoint.partition.eq("standardized_exclusion_repeat")].copy()
    repeated_splits: dict[str, Any] = {}
    for key, method in {
        "matched_structure_only": "A_structure_only",
        "full_solvai": "F_full_solvai",
    }.items():
        subset = repeated.loc[repeated.method.eq(method)]
        values = subset.groupby("repeat", sort=True).absolute_error.mean()
        if values.index.tolist() != [0, 1, 2, 3, 4]:
            raise AssertionError(f"Incomplete repeat series for {method}")
        repeated_splits[key] = {
            "values_kcal_mol": [float(value) for value in values],
            "mean_kcal_mol": float(values.mean()),
            "sd_kcal_mol": float(values.std(ddof=1)),
            "repeats_below_0_20": int(values.lt(0.20).sum()),
        }

    full = blocks["full_solvai"]
    matched = blocks["matched_structure_only"]
    full_abs = full.absolute_error.to_numpy(float)
    matched_abs = matched.absolute_error.to_numpy(float)
    shuffle_primary = shuffle.loc[shuffle.partition.eq("primary")].copy()
    shuffle_mean = (
        shuffle_primary.groupby("molecule_id", as_index=False)
        .agg(y_true=("y_true", "first"), y_pred=("y_pred", "mean"))
        .sort_values("molecule_id")
    )
    shuffled_abs = np.abs(shuffle_mean.y_true - shuffle_mean.y_pred).to_numpy(float)
    uncertainty = {
        "full_solvai_mae": _bootstrap(full_abs, seed=BOOTSTRAP_SEED),
        "full_vs_matched": _bootstrap(full_abs - matched_abs, seed=BOOTSTRAP_SEED + 1),
        "full_vs_pimd8": _bootstrap(full_abs - np.abs(pimd_error), seed=BOOTSTRAP_SEED + 2),
        "full_vs_mean_shuffled": _bootstrap(full_abs - shuffled_abs, seed=BOOTSTRAP_SEED + 3),
    }

    family = (
        full.groupby("functional_group_family", sort=True)
        .absolute_error.agg(["count", "mean"])
        .reset_index()
    )
    chemistry_family = [
        {
            "family": str(row.functional_group_family),
            "n": int(row["count"]),
            "mae_kcal_mol": float(row["mean"]),
        }
        for _, row in family.iterrows()
    ]

    campaign = pd.read_parquet(root / "results/predictions/campaign_oof.parquet")
    headline = pd.read_parquet(root / "results/predictions/headline_oof.parquet")
    hard = pd.read_parquet(root / "results/predictions/hard_holdout_oof.parquet")
    historical = {
        "physics_mixture_comparator": _legacy_method(
            campaign, "Nested structure-only physics mixture"
        ),
        "fixed_exact_connectivity_teachers": _legacy_method(
            headline, "Fixed narrow response + SMD + ConfSolv response"
        ),
        "nested_exact_connectivity_teachers": _legacy_method(
            headline, "Nested narrow SMD teacher selection"
        ),
        "arrow_only_family_holdout": _legacy_method(
            hard, "SMD response direct + ConfSolv response", "family_holdout"
        ),
        "arrow_only_scaffold_holdout": _legacy_method(
            hard, "SMD response direct + ConfSolv response", "scaffold_holdout"
        ),
        "note": (
            "Historical campaign values are retained for provenance and are not "
            "the primary confirmatory evidence."
        ),
    }

    multilambda = pd.read_csv(root / "results/ablations/multilambda_metrics.csv")
    response_metrics = pd.read_csv(root / "results/ablations/multilambda_response_metrics.csv")
    lambda_teachers = pd.read_parquet(root / "results/ablations/pimd2_multilambda_teacher.parquet")
    lambda_counts = {
        code: int(lambda_teachers[f"success__lambda{code}"].fillna(False).sum())
        for code in ("01", "05", "09")
    }
    lambda_counts["complete_three_point"] = int(lambda_teachers.complete_three_point_curve.sum())

    alternative_supervision: dict[str, float] = {}
    for key, filename in {
        "matched_one_seed_base": "smd_plus_confsolv_response_screen.csv",
        "openfe_diagnostics": "openfe_screen.csv",
        "mlff_hierarchy": "mlff_screen.csv",
        "des370k_water_response": "des370k_screen.csv",
    }.items():
        table = pd.read_csv(root / "results/ablations" / filename)
        alternative_supervision[key] = float(table.mae.iloc[0])

    molsolv_metadata = json.loads(
        (root / "data/manifests/molsolv_smd_teacher_metadata.json").read_text()
    )
    confsolv_metadata = json.loads(
        (root / "data/manifests/confsolv_water_teacher_metadata.json").read_text()
    )
    endpoint_metadata = json.loads(
        (root / "data/manifests/endpoint_label_manifest.json").read_text()
    )
    refit_metadata = {
        name: json.loads(
            (
                root / "results/confirmatory/teacher_refits" / name / "training_metadata.json"
            ).read_text()
        )
        for name in ("combisolv_qm", "molsolv_smd", "confsolv")
    }
    data_counts = {
        "molsolv_source_calculations": int(molsolv_metadata["records_total"]),
        "molsolv_exact_connectivity_filtered_structures": 350_391,
        "molsolv_confirmatory_rows": int(refit_metadata["molsolv_smd"]["rows_after_exclusion"]),
        "confsolv_source_water_conformers": int(confsolv_metadata["raw_water_rows"]),
        "confsolv_exact_connectivity_filtered_connectivities": 39_878,
        "confsolv_model_usable_rows": int(refit_metadata["confsolv"]["original_rows"]),
        "confsolv_confirmatory_rows": int(refit_metadata["confsolv"]["rows_after_exclusion"]),
        "combisolv_confirmatory_rows": int(refit_metadata["combisolv_qm"]["rows_after_exclusion"]),
        "endpoint_public_labels": int(endpoint_metadata["endpoint_selection"]["rows"]),
        "arrow85_connectivities_also_in_freesolv": int(benchmark.freesolv_match_id.notna().sum()),
        "standardized_teacher_rows_removed": 56,
    }

    artifact = joblib.load(root / "models/final/head.joblib")
    card = json.loads((root / "models/final/model_card.json").read_text())
    artifact_metrics = {
        "input": "SMILES only",
        "simulation_at_inference": False,
        "descriptor_features": len(artifact["descriptor_columns"]),
        "physics_response_features": len(artifact["teacher_columns"]),
        "total_features": len(artifact["descriptor_columns"]) + len(artifact["teacher_columns"]),
        "ensemble_members": len(artifact["models"]),
        "seeds": [int(seed) for seed in artifact["seeds"]],
        "standardized_exclusion_artifact": bool(card.get("standardized_teacher_exclusions")),
    }

    expected = {
        "matched_structure_only": 0.3033484655954973,
        "full_solvai": 0.20223406721910991,
        "arrow_pimd8": 0.20483647058823526,
        "zero_arrow_full_solvai": 0.25694059186637636,
    }
    for key, value in expected.items():
        observed = float(methods[key]["mae_kcal_mol"])
        if not np.isclose(observed, value, atol=1e-12):
            raise AssertionError(f"Confirmatory metric drift for {key}: {observed} != {value}")
    if artifact_metrics["total_features"] != 2280:
        raise AssertionError("Packaged feature schema drift")
    if not artifact_metrics["standardized_exclusion_artifact"]:
        raise AssertionError("Final artifact is not the standardized-exclusion refit")
    endpoint_audit = next(
        row for row in chemistry_audit["sources"] if row["source"] == "endpoint_experimental"
    )
    if any(
        endpoint_audit[key]
        for key in (
            "full_inchi_key_matches",
            "connectivity_matches",
            "fragment_parent_matches",
            "uncharged_parent_matches",
            "canonical_tautomer_matches",
        )
    ):
        raise AssertionError("Endpoint standardized identity overlap detected")

    metrics: dict[str, Any] = {
        "schema_version": 2,
        "benchmark": {
            "name": "ARROW 85-solute neutral-hydration reference set",
            "molecules": 85,
            "temperature_k": 298.0,
            "freesolv_connectivity_overlap": int(benchmark.freesolv_match_id.notna().sum()),
        },
        "methods": methods,
        "repeated_splits": repeated_splits,
        "global_separation": separation_metrics,
        "paired_confirmatory": paired.to_dict(orient="records"),
        "bootstrap": uncertainty,
        "chemistry_family": chemistry_family,
        "historical_campaign": historical,
        "alternative_supervision": alternative_supervision,
        "multilambda": {
            "teacher_counts": lambda_counts,
            "method_mae_kcal_mol": {
                str(row.method): float(row.mae) for row in multilambda.itertuples(index=False)
            },
            "response_head_metrics": response_metrics.replace({np.nan: None}).to_dict(
                orient="records"
            ),
        },
        "data_counts": data_counts,
        "artifact": artifact_metrics,
        "audit": {
            "assertions_passed": True,
            "standardized_endpoint_matches": 0,
            "standardized_teacher_rows_removed": 56,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "bootstrap_seed": BOOTSTRAP_SEED,
        },
    }

    tidy_rows: list[dict[str, Any]] = []
    for key, value in methods.items():
        for metric in ("mae_kcal_mol", "rmse_kcal_mol"):
            tidy_rows.append(
                {
                    "section": "method",
                    "metric": f"{key}.{metric.removesuffix('_kcal_mol')}",
                    "value": value[metric],
                    "unit": "kcal/mol",
                    "n": value["n"],
                }
            )
    for key, value in repeated_splits.items():
        for metric in ("mean_kcal_mol", "sd_kcal_mol"):
            tidy_rows.append(
                {
                    "section": "repeat",
                    "metric": f"{key}.{metric}",
                    "value": value[metric],
                    "unit": "kcal/mol",
                    "n": 5,
                }
            )
    for regime, regime_values in separation_metrics.items():
        for method, value in regime_values.items():
            tidy_rows.append(
                {
                    "section": "global_separation",
                    "metric": f"{regime}.{method}.mae",
                    "value": value["mae_kcal_mol"],
                    "unit": "kcal/mol",
                    "n": value["n"],
                }
            )
    for key, value in data_counts.items():
        tidy_rows.append(
            {
                "section": "data",
                "metric": key,
                "value": value,
                "unit": "count",
                "n": np.nan,
            }
        )
    return metrics, pd.DataFrame(tidy_rows)


def freeze_markdown(metrics: dict[str, Any]) -> str:
    methods = metrics["methods"]
    repeats = metrics["repeated_splits"]
    data = metrics["data_counts"]
    return f"""# SolvAI paper freeze

This is the canonical quantitative state for the Nature Communications manuscript.
It is generated by `solv_ai.paper_metrics`; hand-edited numerical values are not
authoritative.

## Confirmatory scientific result

Under the preregistered matched endpoint, structure alone gives
**{methods["matched_structure_only"]["mae_kcal_mol"]:.5f} kcal/mol** and the 15 aligned
response priors give **{methods["full_solvai"]["mae_kcal_mol"]:.5f} kcal/mol**. Across
five complete partitions, full SolvAI gives
**{repeats["full_solvai"]["mean_kcal_mol"]:.5f} ±
{repeats["full_solvai"]["sd_kcal_mol"]:.5f} kcal/mol**. The reconstructed
ARROW/PIMD8 reference is **{methods["arrow_pimd8"]["mae_kcal_mol"]:.5f} kcal/mol**.

The conservative headline is therefore: SolvAI reaches the ARROW/PIMD8 accuracy
scale on this reference set while using structure-only inference. Robust sub-0.20
performance and superiority over PIMD8 are not claimed.

## Global separation and transfer

| Quantity | MAE (kcal/mol) |
| --- | ---: |
| Matched structure-only, fixed OOF | {methods["matched_structure_only"]["mae_kcal_mol"]:.5f} |
| Full SolvAI, fixed OOF | {methods["full_solvai"]["mae_kcal_mol"]:.5f} |
| Full SolvAI, five-repeat mean | {repeats["full_solvai"]["mean_kcal_mol"]:.5f} |
| Full SolvAI, no ARROW labels | {methods["zero_arrow_full_solvai"]["mae_kcal_mol"]:.5f} |
| Full SolvAI, global family separation | {metrics["global_separation"]["global_family"]["F_full_solvai"]["mae_kcal_mol"]:.5f} |
| Full SolvAI, global scaffold separation | {metrics["global_separation"]["global_scaffold"]["F_full_solvai"]["mae_kcal_mol"]:.5f} |
| ARROW/PIMD8 | {methods["arrow_pimd8"]["mae_kcal_mol"]:.5f} |

## Data and artifact integrity

- Benchmark molecules: {metrics["benchmark"]["molecules"]}.
- External endpoint labels: {data["endpoint_public_labels"]}.
- Standardized-equivalent teacher rows removed: {data["standardized_teacher_rows_removed"]}.
- Inference features: {metrics["artifact"]["descriptor_features"]} structure +
  {metrics["artifact"]["physics_response_features"]} predicted response priors.
- Simulation at inference: no.
- Every confirmatory prediction table contains complete molecule-level outputs.
- Bootstrap: {BOOTSTRAP_REPLICATES:,} paired molecule resamples, seed
  {BOOTSTRAP_SEED}.

Historical 0.238606 and 0.197047 campaign values remain in the machine-readable
`historical_campaign` section for provenance but are not the matched confirmatory
headline.
"""


def write_paper_metrics(root: Path = ROOT) -> dict[str, Any]:
    metrics, tidy = compute_paper_metrics(root)
    serialized = json.dumps(metrics, indent=2) + "\n"
    (root / "results/paper_metrics.json").write_text(serialized)
    # Preserve the requested deliverable name without maintaining a second,
    # potentially stale scientific record.
    (root / "results/final_metrics.json").write_text(serialized)
    tidy.to_csv(root / "results/paper_metrics.csv", index=False)
    (root / "reports/PAPER_FREEZE.md").write_text(freeze_markdown(metrics))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    metrics = write_paper_metrics(args.root.resolve())
    print(json.dumps(metrics["methods"], indent=2))


if __name__ == "__main__":
    main()
