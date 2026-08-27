"""Write the final structure-only physics-distillation metrics and report."""

from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arrow_distill.data import ROOT


def main() -> None:
    confirmation = pd.read_csv(
        ROOT / "results/nested_smd_confsolv_structure_confirmation.csv"
    ).set_index("method")
    hard = pd.read_csv(ROOT / "results/smd_confsolv_response_structure_hard.csv")
    repeated = pd.read_csv(
        ROOT / "results/repeated_nested_smd_confsolv_confirmation_summary.csv"
    ).set_index("method")
    repeated_detail = json.loads(
        (ROOT / "results/repeated_structure_only_confirmation.json").read_text()
    )
    master = pd.read_parquet(ROOT / "data/processed/arrow_solvation_master.parquet")
    hydration = master[master.solvent.eq("water")]
    raw_pimd8 = float((hydration.delta_g_pimd8 - hydration.delta_g_exp).abs().mean())
    fixed_name = "Fixed narrow response + SMD + ConfSolv response"
    nested_name = "Nested narrow SMD teacher selection"
    repeated_fixed_name = "narrow response + SMD + ConfSolv response"
    fixed = float(confirmation.loc[fixed_name, "mae"])
    nested = float(confirmation.loc[nested_name, "mae"])
    family = float(hard.loc[hard.regime.eq("family_holdout"), "mae"].iloc[0])
    scaffold = float(hard.loc[hard.regime.eq("scaffold_holdout"), "mae"].iloc[0])
    repeat_fixed = float(repeated.loc[repeated_fixed_name, "mean_mae"])
    repeat_fixed_std = float(repeated.loc[repeated_fixed_name, "std_mae"])
    repeat_nested = float(repeated.loc["Nested selection", "mean_mae"])

    rows = [
        {
            "method": "ARROW/PIMD8",
            "simulation_at_inference": True,
            "random_oof_mae": raw_pimd8,
            "repeat_mean_mae": np.nan,
            "family_mae": np.nan,
            "scaffold_mae": np.nan,
            "status": "reference",
        },
        {
            "method": "PIMD8 + nested ML residual",
            "simulation_at_inference": True,
            "random_oof_mae": 0.18682171894402116,
            "repeat_mean_mae": np.nan,
            "family_mae": 0.19644319866557736,
            "scaffold_mae": 0.20055713815126047,
            "status": "ceiling only",
        },
        {
            "method": "Previous structure-only baseline",
            "simulation_at_inference": False,
            "random_oof_mae": 0.23860611898039194,
            "repeat_mean_mae": np.nan,
            "family_mae": 0.31698158379458463,
            "scaffold_mae": 0.3289238985826331,
            "status": "superseded",
        },
        {
            "method": "Structure + narrow response",
            "simulation_at_inference": False,
            "random_oof_mae": float(confirmation.loc["Fixed narrow response without SMD", "mae"]),
            "repeat_mean_mae": float(repeated.loc["narrow response without SMD", "mean_mae"]),
            "family_mae": np.nan,
            "scaffold_mae": np.nan,
            "status": "ablation",
        },
        {
            "method": "Structure + SMD teacher",
            "simulation_at_inference": False,
            "random_oof_mae": float(confirmation.loc["Fixed narrow response + SMD water", "mae"]),
            "repeat_mean_mae": float(repeated.loc["narrow response + SMD water", "mean_mae"]),
            "family_mae": np.nan,
            "scaffold_mae": np.nan,
            "status": "ablation",
        },
        {
            "method": "Fixed SMD + ConfSolv student",
            "simulation_at_inference": False,
            "random_oof_mae": fixed,
            "repeat_mean_mae": repeat_fixed,
            "family_mae": family,
            "scaffold_mae": scaffold,
            "status": "packaged candidate",
        },
        {
            "method": "Nested SMD/ConfSolv selection",
            "simulation_at_inference": False,
            "random_oof_mae": nested,
            "repeat_mean_mae": repeat_nested,
            "family_mae": np.nan,
            "scaffold_mae": np.nan,
            "status": "selection-adjusted headline",
        },
    ]
    screens = (
        ("smd_confsolv_openfe_diagnostics_screen.csv", "OpenFE diagnostics"),
        ("smd_confsolv_mlff_hierarchy_screen.csv", "MLFF HFE hierarchy"),
        ("smd_confsolv_des_water_screen.csv", "DES370K water response"),
    )
    for filename, name in screens:
        path = ROOT / "results" / filename
        if path.is_file():
            rows.append(
                {
                    "method": name,
                    "simulation_at_inference": False,
                    "random_oof_mae": float(pd.read_csv(path).mae.iloc[0]),
                    "repeat_mean_mae": np.nan,
                    "family_mae": np.nan,
                    "scaffold_mae": np.nan,
                    "status": "one-seed new-information screen",
                }
            )
    lambda_path = ROOT / "results/multilambda_physics_distillation_oof.csv"
    lambda_metrics = pd.read_csv(lambda_path) if lambda_path.is_file() else pd.DataFrame()
    for row in lambda_metrics.itertuples(index=False):
        rows.append(
            {
                "method": row.method,
                "simulation_at_inference": False,
                "random_oof_mae": float(row.mae),
                "repeat_mean_mae": np.nan,
                "family_mae": np.nan,
                "scaffold_mae": np.nan,
                "status": "fixed multi-lambda ablation",
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(ROOT / "results/physics_distillation_model_comparison.csv", index=False)

    curve_path = ROOT / "results/multilambda_response_student_metrics.csv"
    curve_metrics = pd.read_csv(curve_path) if curve_path.is_file() else pd.DataFrame()
    teacher_path = ROOT / "results/pimd2_multilambda_teacher.parquet"
    teachers = pd.read_parquet(teacher_path) if teacher_path.is_file() else pd.DataFrame()
    complete_curves = int(teachers.complete_three_point_curve.sum()) if len(teachers) else 0
    teacher_counts = (
        {
            label: int(teachers[f"success__{label}"].fillna(False).sum())
            for label in ("lambda01", "lambda05", "lambda09")
        }
        if len(teachers)
        else {}
    )
    repeated_stats = repeated_detail["methods"][repeated_fixed_name]
    final_metrics = {
        "best_fixed_structure_only_random_oof_mae": fixed,
        "selection_adjusted_nested_random_oof_mae": nested,
        "fixed_structure_only_repeat_mean_mae": repeat_fixed,
        "fixed_structure_only_repeat_std_mae": repeat_fixed_std,
        "nested_repeat_mean_mae": repeat_nested,
        "fixed_repeat_molecule_bootstrap_95_ci": repeated_stats["molecule_cluster_bootstrap_95_ci"],
        "fixed_repeat_bootstrap_probability_below_0_20": repeated_stats[
            "bootstrap_probability_mae_below_0_20"
        ],
        "family_holdout_mae": family,
        "scaffold_holdout_mae": scaffold,
        "raw_pimd8_mae": raw_pimd8,
        "pimd8_nested_residual_mae": 0.18682171894402116,
        "pimd8_labels_used_in_candidate_distillation": 85,
        "pimd8_labels_used_by_selected_artifact": 0,
        "complete_three_point_pimd2_teachers": complete_curves,
        "pimd2_teacher_success_counts": teacher_counts,
        "single_split_below_0_20": bool(nested < 0.20),
        "repeat_mean_below_0_20": bool(repeat_fixed < 0.20),
        "family_holdout_below_0_20": bool(family < 0.20),
        "inference_simulation": False,
        "new_pimd8_simulations": 0,
    }
    if len(lambda_metrics):
        final_metrics["multilambda_ablation"] = dict(
            zip(lambda_metrics.method, lambda_metrics.mae, strict=True)
        )
        short_names = {
            "Multi-lambda physics distillation A: structure/response baseline": "Matched baseline",
            "Multi-lambda physics distillation B2: +distilled PIMD2 lambda response": "+ predicted lambda response",
            "Multi-lambda physics distillation B1: +distilled classical-NQE-PIMD hierarchy": "+ classical/NQE/PIMD hierarchy",
            "Multi-lambda physics distillation B: +full distilled physics hierarchy": "+ both physics blocks",
            "Multi-lambda physics distillation C: integrated predicted dH/dlambda (affine calibration)": "Integrated predicted curve",
        }
        plot_frame = (
            lambda_metrics.loc[lambda_metrics.mae.lt(0.5)]
            .assign(label=lambda frame: frame.method.map(short_names).fillna(frame.method))
            .sort_values("mae", ascending=True)
        )
        figure, axis = plt.subplots(figsize=(8.5, 4.6))
        colors = ["#2a9d8f" if value < 0.20 else "#e76f51" for value in plot_frame.mae]
        axis.barh(plot_frame.label, plot_frame.mae, color=colors)
        axis.axvline(0.20, color="black", linestyle="--", linewidth=1.2, label="0.20 target")
        axis.set_xlabel("Strict OOF MAE (kcal/mol)")
        axis.set_title("Structure-only multi-lambda distillation ablation")
        axis.set_xlim(0.185, 0.225)
        axis.legend(loc="lower right")
        figure.tight_layout()
        figure.savefig(ROOT / "figures/structure_only_multilambda_ablation.png", dpi=220)
        plt.close(figure)
    (ROOT / "results/physics_distillation_final_metrics.json").write_text(
        json.dumps(final_metrics, indent=2) + "\n"
    )

    lambda_lines = (
        "\n".join(
            f"- {row.method}: **{row.mae:.5f} kcal/mol**"
            for row in lambda_metrics.itertuples(index=False)
        )
        if len(lambda_metrics)
        else "- Final multi-lambda run pending."
    )
    response_summary = (
        f"Nonconstant response-head MAE ranges from "
        f"{curve_metrics.loc[curve_metrics.target_std.gt(0), 'mae'].min():.2f} to "
        f"{curve_metrics.loc[curve_metrics.target_std.gt(0), 'mae'].max():.2f} "
        "kcal/mol."
        if len(curve_metrics)
        else "Response-head diagnostics pending."
    )
    report = f"""# Structure-only ARROW/PIMD physics distillation: final report

## Decision

The structure-only model crosses 0.20 on one leakage-safe five-fold partition:
the fixed SMD+ConfSolv student reaches **{fixed:.5f} kcal/mol**, and nested
feature-block selection reaches **{nested:.5f}**. The deployed artifact consumes
SMILES only and runs no simulation.

The crossing is **not robustly confirmed**. Across five independent outer-CV
partitions, the fixed student averages **{repeat_fixed:.5f} ±
{repeat_fixed_std:.5f}**, nested selection averages **{repeat_nested:.5f}**, and
the fixed student's molecule-cluster bootstrap 95% interval is
**[{repeated_stats["molecule_cluster_bootstrap_95_ci"][0]:.5f},
{repeated_stats["molecule_cluster_bootstrap_95_ci"][1]:.5f}]**. Family and scaffold
holdouts remain **{family:.5f}** and **{scaffold:.5f}**. The scientifically honest
answer is therefore: a sub-0.20 point estimate exists, but robust sub-0.20
generalization has not been demonstrated.

## Headline comparison

| Method | Simulation at inference? | Fixed random OOF | Five-repeat mean | Family | Scaffold |
|---|---:|---:|---:|---:|---:|
| ARROW/PIMD8 | Yes | {raw_pimd8:.5f} | — | — | — |
| PIMD8 + nested residual | Yes | 0.18682 | — | 0.19644 | 0.20056 |
| Previous structure-only baseline | No | 0.23861 | — | 0.31698 | 0.32892 |
| Structure + narrow response | No | {confirmation.loc["Fixed narrow response without SMD", "mae"]:.5f} | {repeated.loc["narrow response without SMD", "mean_mae"]:.5f} | — | — |
| + MolSolv SMD teacher | No | {confirmation.loc["Fixed narrow response + SMD water", "mae"]:.5f} | {repeated.loc["narrow response + SMD water", "mean_mae"]:.5f} | — | — |
| + ConfSolv response (fixed) | No | **{fixed:.5f}** | **{repeat_fixed:.5f}** | {family:.5f} | {scaffold:.5f} |
| Nested SMD/ConfSolv selection | No | **{nested:.5f}** | {repeat_nested:.5f} | — | — |

## What transferred

MolSolv contributes 350,391 benchmark-disjoint SMD(water) structures; adding its
water-response teacher improves the fixed model from 0.21362 to 0.20161. ConfSolv
adds 39,878 benchmark-disjoint water structures with conformer-response moments and
improves the fixed result further to 0.19705. This is the strongest direct evidence
that physics-rich supervision transfers into a structure-only student.

Broader physical blocks did not help this already strong representation. Predicted
OpenFE diagnostics, the MLFF/force-field hierarchy, and DES370K water/SAPT response
score 0.20575, 0.20109, and 0.20993 in matched one-seed screens versus 0.19191 for
their unchanged base. They were stopped rather than combined or tuned.

## Lambda-response experiment

{lambda_lines}

{response_summary} All measured PIMD2 values are training-only labels. For every
outer-test molecule, both its response curve and classical/PIMD hierarchy labels
are withheld; the final prediction uses only its structure.

## Leakage and deployment

The external sources used by the final artifact have zero connectivity overlap
with all 85 benchmark molecules. The packaged feature audit finds zero forbidden
test-time fields. End-to-end inference reconstructs all 15 learned physics features
from SMILES within 1.9e-6 of their cached values. The artifact is under
`models/final_structure_only/`; the CLI is `scripts/predict_structure_only.py`.

## Remaining blocker and exact next data

The limiting step is not the final regressor. It is structure-to-response transfer
from sparse, chemically narrow high-fidelity labels. The largest errors are amides,
aromatics, ethers, acids, and alkanes, while current independent classical/PIMD or
lambda-response sets have low nearest-neighbor similarity and only tens of examples
per family. The next acquisition should be a protocol-matched, benchmark-disjoint
set containing full lambda-resolved dH/dlambda plus electrostatic, polarization,
dispersion/repulsion, and classical/PIMD8 pairs for at least 50 diverse molecules
in each of those five families (roughly 250-500 molecules total). That directly
targets the observed response-surrogate error; more generic SMILES tuning does not.

## Required conclusions

1. Best zero-simulation experimental MAE: {fixed:.5f} fixed OOF; {nested:.5f} selection-adjusted nested OOF.
2. Best family-held-out MAE: {family:.5f}.
3. Raw ARROW/PIMD8 MAE on comparable data: {raw_pimd8:.5f}.
4. PIMD8 + ML residual MAE: 0.18682 (not eligible for deployment).
5. PIMD-distilled MAE: 0.29169 for the earlier PIMD-only student.
6. Final ensemble MAE: {fixed:.5f} on the original split; {repeat_fixed:.5f} five-repeat mean.
7. Number of PIMD8 labels used: 85 in candidate distillation experiments; 0 in the selected artifact because that block worsened OOF validation.
8. Number of NEW PIMD8 simulations performed: 0.
9. Estimated inference speedup vs PIMD8: not used as a claim; inference is structure-model evaluation only.
10. Does the project beat 0.20 legitimately? NO — one strict split crosses, but repeated/hard validation does not.
11. Does it beat 0.20 without simulation at inference? NO robustly; YES only for the original strict nested point estimate.
12. Most important model component: benchmark-disjoint SMD(water) plus ConfSolv conformer-response distillation.
13. Main failure family: amides.
14. Best next experiment if given one additional day: acquire protocol-matched full lambda/component and paired classical/PIMD labels for the five high-error families above.
"""
    (ROOT / "reports/PHYSICS_DISTILLATION_FINAL_REPORT.md").write_text(report)
    print(json.dumps(final_metrics, indent=2))


if __name__ == "__main__":
    main()
