"""Derive collaborator follow-up diagnostics from frozen SolvAI outputs.

This script does not fit or select a model. It summarizes the standardized-
equivalence confirmatory predictions, predeclared source-block controls and frozen
teacher validation records.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "followup_audit"
WORKSPACE = ROOT.parents[1]
BOOTSTRAP_SEED = 20260828


def paired_interval(
    truth: np.ndarray, candidate: np.ndarray, baseline: np.ndarray, draws: int = 100_000
) -> tuple[float, float, float, float]:
    delta = np.abs(truth - candidate) - np.abs(truth - baseline)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, len(delta), size=(draws, len(delta)))
    sampled = delta[indices].mean(axis=1)
    low, high = np.quantile(sampled, [0.025, 0.975])
    return float(delta.mean()), float(low), float(high), float(np.mean(delta < 0))


def source_blocks() -> pd.DataFrame:
    metrics = pd.read_csv(ROOT / "results/confirmatory/standardized_exclusion_endpoint_metrics.csv")
    predictions = pd.read_parquet(
        ROOT / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    primary_name = "standardized_exclusion_primary"
    primary = predictions.loc[predictions.partition.eq(primary_name)]
    base = primary.loc[primary.method.eq("A_structure_only")].sort_values("molecule_id")
    base_mae = float(np.mean(base.absolute_error))
    descriptions = {
        "A_structure_only": ("Structure only", "none", "matched causal baseline"),
        "B_empirical_residual": (
            "Empirical/residual-corrected",
            "Abraham E/S/A/B/L + corrected OpenFF + corrected GBn2",
            "strongest isolated block",
        ),
        "C_computation_core": (
            "Computation-derived core",
            "CombiSolv-QM + raw OpenFF + raw GBn2",
            "lower point estimate alone; paired interval crosses zero",
        ),
        "D_smd_water": (
            "SMD(water)",
            "MolSolv SMD(water)",
            "lower point estimate alone; paired interval crosses zero",
        ),
        "E_confsolv": (
            "ConfSolv",
            "six conformational/water-response summaries",
            "neutral alone; complementary in stack",
        ),
        "G_narrow_reference": (
            "Narrow cumulative stack",
            "CombiSolv-QM + Abraham + corrected OpenFF/GBn2",
            "secondary cumulative comparison",
        ),
        "H_narrow_smd_reference": (
            "Narrow + SMD",
            "narrow stack + MolSolv SMD(water)",
            "secondary cumulative comparison",
        ),
        "F_full_solvai": (
            "Full SolvAI",
            "all 15 response priors",
            "best predeclared stack; supports complementarity",
        ),
    }
    rows: list[dict[str, object]] = []
    for method, (label, contents, interpretation) in descriptions.items():
        selected = primary.loc[primary.method.eq(method)].sort_values("molecule_id")
        if len(selected) != 85:
            raise AssertionError(f"Expected 85 primary predictions for {method}")
        if method == "A_structure_only":
            delta, low, high, fraction = 0.0, np.nan, np.nan, np.nan
        else:
            delta, low, high, fraction = paired_interval(
                selected.y_true.to_numpy(),
                selected.y_pred.to_numpy(),
                base.y_pred.to_numpy(),
            )
        metric = metrics.loc[metrics.partition.eq(primary_name) & metrics.method.eq(method)].iloc[0]
        rows.append(
            {
                "method": method,
                "label": label,
                "contents": contents,
                "mae_kcal_mol": metric.mae,
                "rmse_kcal_mol": metric.rmse,
                "mae_change_vs_structure": delta,
                "paired_ci_low_95": low,
                "paired_ci_high_95": high,
                "fraction_molecules_improved": fraction,
                "interpretation": interpretation,
            }
        )
    result = pd.DataFrame(rows)
    if not np.isclose(
        result.loc[result.method.eq("A_structure_only"), "mae_kcal_mol"].iat[0], base_mae
    ):
        raise AssertionError("Structure-only metric mismatch")
    return result


def teacher_fidelity() -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def append(source: str, target: str, n: int, mae: float, units: str, protocol: str) -> None:
        rows.append(
            {
                "source": source,
                "target": target,
                "n_evaluation": n,
                "mae": mae,
                "units": units,
                "evaluation_protocol": protocol,
            }
        )

    for directory, target, source in (
        ("combisolv_qm", "target_qm_water", "CombiSolv-QM"),
        ("molsolv_smd", "target_smd_water", "MolSolv SMD(water)"),
    ):
        base = ROOT / "results/confirmatory/teacher_refits" / directory
        truth = pd.read_csv(base / "input/test.csv")
        prediction = pd.read_csv(base / "model/model_0/test_predictions.csv")
        if not truth.smiles.equals(prediction.smiles):
            raise AssertionError(f"Teacher test ordering mismatch for {source}")
        mae = float(np.mean(np.abs(truth[target] - prediction[target])))
        append(source, target, len(truth), mae, "kcal mol-1", "fixed held-out source test")

    abraham = json.loads((ROOT / "data/manifests/soluteml_hydration_metadata.json").read_text())
    for target, values in abraham["abraham_targets"].items():
        append(
            "SoluteML Abraham",
            target,
            int(values["n_source"]),
            float(values["oof_mae"]),
            "Abraham scale",
            "source OOF",
        )

    for filename, source, target_keys in (
        (
            "openff_alchemical_teacher_metadata.json",
            "OpenFF explicit-water",
            ("openff23_dg", "openff23_exp_residual"),
        ),
        (
            "implicit_solvent_teacher_metadata.json",
            "GBn2 implicit-solvent",
            ("gbn2_alchemical_dg", "gbn2_exp_residual"),
        ),
    ):
        metadata = json.loads((ROOT / "data/manifests" / filename).read_text())
        for target in target_keys:
            values = metadata["targets"][target]
            append(
                source,
                target,
                int(values["n_source"]),
                float(values["oof_mae"]),
                "kcal mol-1",
                "source OOF",
            )

    confsolv = json.loads(
        (ROOT / "results/confirmatory/teacher_refits/confsolv/training_metadata.json").read_text()
    )
    selected = {
        "confsolv_gas_conformer_correction",
        "confsolv_solution_conformer_correction",
        "confsolv_hydration_conformer_correction",
        "confsolv_water_gsolv_std",
        "confsolv_water_response_mean",
        "confsolv_water_response_std",
    }
    for target, values in confsolv["validation"].items():
        if target in selected:
            append(
                "ConfSolv H2O",
                target,
                int(values["n_validation"]),
                float(values["mae_kcal_mol"]),
                "kcal mol-1",
                "fixed source validation",
            )
    return pd.DataFrame(rows)


def hard_cases() -> tuple[pd.DataFrame, pd.DataFrame]:
    predictions = pd.read_parquet(
        ROOT / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    primary = predictions.loc[
        predictions.partition.eq("standardized_exclusion_primary")
        & predictions.method.isin(["A_structure_only", "F_full_solvai"])
    ]
    wide = primary.pivot(
        index="molecule_id",
        columns="method",
        values=["y_true", "y_pred", "absolute_error"],
    )
    metadata = primary.drop_duplicates("molecule_id").set_index("molecule_id")[
        [
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
        ]
    ]
    benchmark = pd.read_parquet(ROOT / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")].drop_duplicates("molecule_id")
    benchmark = benchmark.set_index("molecule_id")
    columns = [
        "delta_g_pimd8",
        "delta_g_classical_arrow",
        "heavy_atom_count",
        "molecular_weight",
        "hbd",
        "hba",
        "tpsa",
        "logp",
        "rotatable_bonds",
    ]
    result = metadata.join(benchmark[columns], how="left")
    result["delta_g_exp"] = wide["y_true"]["F_full_solvai"]
    result["delta_g_solvai"] = wide["y_pred"]["F_full_solvai"]
    result["solvai_absolute_error"] = wide["absolute_error"]["F_full_solvai"]
    result["delta_g_structure_only"] = wide["y_pred"]["A_structure_only"]
    result["structure_only_absolute_error"] = wide["absolute_error"]["A_structure_only"]
    result["absolute_error_improvement"] = (
        result.structure_only_absolute_error - result.solvai_absolute_error
    )
    result["pimd8_absolute_error"] = np.abs(result.delta_g_pimd8 - result.delta_g_exp)
    result["solvai_signed_error"] = result.delta_g_solvai - result.delta_g_exp
    result["pimd8_signed_error"] = result.delta_g_pimd8 - result.delta_g_exp
    result["solvai_pimd8_error_same_direction"] = np.sign(result.solvai_signed_error) == np.sign(
        result.pimd8_signed_error
    )
    result["pimd8_more_accurate"] = result.pimd8_absolute_error < result.solvai_absolute_error
    result["solvai_large_pimd8_small"] = result.solvai_absolute_error.ge(
        0.5
    ) & result.pimd8_absolute_error.lt(0.5)
    aromatic_rings = []
    total_rings = []
    for smiles in result.canonical_smiles:
        molecule = Chem.MolFromSmiles(smiles)
        aromatic_rings.append(rdMolDescriptors.CalcNumAromaticRings(molecule))
        total_rings.append(rdMolDescriptors.CalcNumRings(molecule))
    result["aromatic_ring_count"] = aromatic_rings
    result["ring_count"] = total_rings
    result["larger_aromatic_posthoc"] = result.aromatic_ring_count.gt(
        0
    ) & result.heavy_atom_count.ge(8)

    repeated = predictions.loc[
        predictions.partition.eq("standardized_exclusion_repeat")
        & predictions.method.eq("F_full_solvai")
    ]
    repeated_summary = repeated.groupby("molecule_id").agg(
        repeat_mean_prediction=("y_pred", "mean"),
        repeat_prediction_sd=("y_pred", "std"),
        repeat_mean_absolute_error=("absolute_error", "mean"),
        repeat_max_absolute_error=("absolute_error", "max"),
    )
    result = result.join(repeated_summary, how="left")
    result = result.sort_values("solvai_absolute_error", ascending=False).reset_index()
    result.insert(0, "solvai_error_rank", np.arange(1, len(result) + 1))

    groups: list[dict[str, object]] = []
    masks = {
        "all": np.ones(len(result), dtype=bool),
        "amides": result.functional_group_family.eq("Amides").to_numpy(),
        "all_aromatic_structures": result.aromatic_ring_count.gt(0).to_numpy(),
        "larger_aromatics_heavy_atoms_ge_8_posthoc": result.larger_aromatic_posthoc.to_numpy(),
        "all_heavy_atoms_ge_8_posthoc": result.heavy_atom_count.ge(8).to_numpy(),
    }
    for label, mask in masks.items():
        subset = result.loc[mask]
        groups.append(
            {
                "diagnostic_group": label,
                "n": len(subset),
                "solvai_mae": subset.solvai_absolute_error.mean(),
                "structure_only_mae": subset.structure_only_absolute_error.mean(),
                "pimd8_mae": subset.pimd8_absolute_error.mean(),
                "mean_heavy_atoms": subset.heavy_atom_count.mean(),
                "posthoc": label != "all",
            }
        )
    for family, subset in result.groupby("functional_group_family", sort=False):
        groups.append(
            {
                "diagnostic_group": f"family::{family}",
                "n": len(subset),
                "solvai_mae": subset.solvai_absolute_error.mean(),
                "structure_only_mae": subset.structure_only_absolute_error.mean(),
                "pimd8_mae": subset.pimd8_absolute_error.mean(),
                "mean_heavy_atoms": subset.heavy_atom_count.mean(),
                "posthoc": True,
            }
        )
    return result, pd.DataFrame(groups)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_blocks().to_csv(OUT / "source_block_summary.csv", index=False)
    teacher_fidelity().to_csv(OUT / "teacher_fidelity_summary.csv", index=False)
    hard, groups = hard_cases()
    hard.to_csv(OUT / "hard_case_audit.csv", index=False)
    groups.to_csv(OUT / "hard_case_group_summary.csv", index=False)
    metadata = {
        "analysis_type": "descriptive analysis of frozen outputs; no model fitting or selection",
        "primary_predictions": "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "outputs": [
            "source_block_summary.csv",
            "teacher_fidelity_summary.csv",
            "hard_case_audit.csv",
            "hard_case_group_summary.csv",
        ],
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
