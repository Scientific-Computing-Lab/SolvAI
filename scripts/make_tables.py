#!/usr/bin/env python3
"""Render manuscript tables and TeX macros from the frozen metric file."""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paper/tables"


def main() -> None:
    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    methods = metrics["methods"]
    repeats = metrics["repeated_splits"]
    required = {
        "previous_structure_only": 0.23860611898039202,
        "smd_confsolv_fixed": 0.19704747409482312,
        "nested_selection": 0.19930603930646335,
        "arrow_pimd8": 0.20483647058823526,
    }
    for key, expected in required.items():
        observed = methods[key]["mae_kcal_mol"]
        if abs(observed - expected) > 1e-12:
            raise AssertionError(f"Refusing to render stale metric {key}: {observed}")

    OUT.mkdir(parents=True, exist_ok=True)
    macros = {
        "BenchmarkN": str(metrics["benchmark"]["molecules"]),
        "ClassicalMAE": f"{methods['classical_arrow']['mae_kcal_mol']:.3f}",
        "PIMDMAE": f"{methods['arrow_pimd8']['mae_kcal_mol']:.3f}",
        "PreviousMAE": f"{methods['previous_structure_only']['mae_kcal_mol']:.3f}",
        "NarrowMAE": f"{methods['narrow_response']['mae_kcal_mol']:.3f}",
        "SMDMAE": f"{methods['smd_water']['mae_kcal_mol']:.3f}",
        "SolvAIMAE": f"{methods['smd_confsolv_fixed']['mae_kcal_mol']:.3f}",
        "NestedMAE": f"{methods['nested_selection']['mae_kcal_mol']:.3f}",
        "RepeatMean": f"{repeats['fixed']['mean_kcal_mol']:.3f}",
        "RepeatSD": f"{repeats['fixed']['sd_kcal_mol']:.3f}",
        "NestedRepeatMean": f"{repeats['nested']['mean_kcal_mol']:.3f}",
        "NestedRepeatSD": f"{repeats['nested']['sd_kcal_mol']:.3f}",
        "FamilyMAE": f"{methods['family_holdout']['mae_kcal_mol']:.3f}",
        "ScaffoldMAE": f"{methods['scaffold_holdout']['mae_kcal_mol']:.3f}",
        "MolSolvN": f"{metrics['data_counts']['molsolv_training_structures']:,}",
        "ConfSolvN": f"{metrics['data_counts']['confsolv_training_connectivities']:,}",
    }
    (OUT / "metrics_macros.tex").write_text(
        "".join(f"\\newcommand{{\\{name}}}{{{value}}}\n" for name, value in macros.items())
    )

    rows = [
        ("Classical ARROW", "yes", "classical_arrow", "--", "--"),
        ("ARROW/PIMD8", "yes", "arrow_pimd8", "--", "--"),
        ("Previous structure-only", "no", "previous_structure_only", "--", "--"),
        ("Narrow response", "no", "narrow_response", "--", "--"),
        ("+ SMD(water)", "no", "smd_water", "--", "--"),
        (
            "SolvAI: + ConfSolv response",
            "no",
            "smd_confsolv_fixed",
            "family_holdout",
            "scaffold_holdout",
        ),
        ("Nested feature selection", "no", "nested_selection", "--", "--"),
    ]
    lines = [
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Method & Simulation at inference & Random OOF & Family holdout & Scaffold holdout \\",
        r"\midrule",
    ]
    csv_rows = []
    for label, simulation, metric_key, family_key, scaffold_key in rows:
        random = methods[metric_key]["mae_kcal_mol"]
        family = "--" if family_key == "--" else f"{methods[family_key]['mae_kcal_mol']:.3f}"
        scaffold = "--" if scaffold_key == "--" else f"{methods[scaffold_key]['mae_kcal_mol']:.3f}"
        lines.append(f"{label} & {simulation} & {random:.3f} & {family} & {scaffold} \\\\")
        csv_rows.append(
            {
                "method": label,
                "simulation_at_inference": simulation,
                "random_oof_mae": random,
                "family_holdout_mae": family,
                "scaffold_holdout_mae": scaffold,
            }
        )
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    (OUT / "main_comparison.tex").write_text("\n".join(lines) + "\n")
    with (OUT / "main_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_rows[0])
        writer.writeheader()
        writer.writerows(csv_rows)


if __name__ == "__main__":
    main()
