#!/usr/bin/env python3
"""Render publication tables and Supplementary Data from frozen artifacts."""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "paper/tables"
SI_TABLES = ROOT / "paper/supplementary/tables"
SUPP_DATA = ROOT / "paper/supplementary_data"
ED = ROOT / "paper/extended_data"


def latex_escape(value: object) -> str:
    text = str(value)
    for source, target in (("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#")):
        text = text.replace(source, target)
    return text


def latex_table(frame: pd.DataFrame, path: Path, widths: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = widths or ("l" * len(frame.columns))
    lines = [rf"\begin{{tabular}}{{{columns}}}", r"\toprule"]
    lines.append(" & ".join(latex_escape(c) for c in frame.columns) + r" \\")
    lines.append(r"\midrule")
    for row in frame.itertuples(index=False, name=None):
        lines.append(" & ".join(latex_escape(value) for value in row) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    path.write_text("\n".join(lines) + "\n")


def workbook(path: Path, sheets: dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        frozen_time = datetime(2026, 8, 27, tzinfo=UTC)
        writer.book.properties.created = frozen_time
        writer.book.properties.modified = frozen_time
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name[:31], index=False)
            ws = writer.book[name[:31]]
            for cell in ws[1]:
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="1769AA")
                cell.alignment = Alignment(wrap_text=True)
            ws.freeze_panes = "A2"
            for column in ws.columns:
                values = [str(c.value or "") for c in column[:1000]]
                ws.column_dimensions[column[0].column_letter].width = min(
                    60, max(10, max(map(len, values)) + 2)
                )
    # OOXML is a ZIP container. Normalize both document properties and member
    # timestamps so Supplementary Data files reproduce byte-for-byte in CI.
    normalized = path.with_name(f".{path.name}.normalized")
    with ZipFile(path) as source, ZipFile(
        normalized, "w", compression=ZIP_DEFLATED, compresslevel=9
    ) as target:
        for name in sorted(source.namelist()):
            payload = source.read(name)
            if name == "docProps/core.xml":
                text = payload.decode("utf-8")
                text = re.sub(
                    r"(<dcterms:(?:created|modified)[^>]*>).*?(</dcterms:(?:created|modified)>)",
                    r"\g<1>2026-08-27T00:00:00Z\g<2>",
                    text,
                )
                payload = text.encode("utf-8")
            source_info = source.getinfo(name)
            info = ZipInfo(name, date_time=(2026, 8, 27, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.external_attr = source_info.external_attr
            info.internal_attr = source_info.internal_attr
            info.create_system = source_info.create_system
            target.writestr(info, payload, compress_type=ZIP_DEFLATED, compresslevel=9)
    normalized.replace(path)


def response_priors() -> pd.DataFrame:
    rows = [
        (
            "combisolv_qm_teacher",
            "COSMOtherm hydration response",
            "CombiSolv-QM water",
            "kcal mol-1",
            "CHEMELEON-initialized D-MPNN",
            3963,
            "direct surrogate",
        ),
        (
            "abraham_e_teacher",
            "excess molar refraction axis E",
            "SoluteML Abraham",
            "Abraham scale",
            "ExtraTrees",
            8098,
            "direct surrogate",
        ),
        (
            "abraham_s_teacher",
            "dipolarity/polarizability axis S",
            "SoluteML Abraham",
            "Abraham scale",
            "ExtraTrees",
            8098,
            "direct surrogate",
        ),
        (
            "abraham_a_teacher",
            "hydrogen-bond acidity axis A",
            "SoluteML Abraham",
            "Abraham scale",
            "ExtraTrees",
            8098,
            "direct surrogate",
        ),
        (
            "abraham_b_teacher",
            "hydrogen-bond basicity axis B",
            "SoluteML Abraham",
            "Abraham scale",
            "ExtraTrees",
            8098,
            "direct surrogate",
        ),
        (
            "abraham_l_teacher",
            "hexadecane-air partition axis L",
            "SoluteML Abraham",
            "Abraham scale",
            "ExtraTrees",
            8098,
            "direct surrogate",
        ),
        (
            "openff_corrected_teacher",
            "explicit-water alchemical hydration response",
            "OpenFF 2.3.0 ASFE",
            "kcal mol-1",
            "ExtraTrees",
            520,
            "predicted calculation + predicted experimental residual",
        ),
        (
            "gbn2_corrected_teacher",
            "implicit-solvent hydration response",
            "GBn2 / GNNImplicitSolvent",
            "kcal mol-1",
            "ExtraTrees",
            550,
            "predicted GBn2 value + predicted experimental residual",
        ),
        (
            "molsolv_smd_teacher",
            "SMD(water) solvation response",
            "MolSolv",
            "kcal mol-1",
            "CHEMELEON-initialized D-MPNN",
            350391,
            "direct surrogate",
        ),
        (
            "confsolv_gas_conformer_correction_teacher",
            "gas conformer free-energy correction",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
        (
            "confsolv_solution_conformer_correction_teacher",
            "solution conformer free-energy correction",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
        (
            "confsolv_hydration_conformer_correction_teacher",
            "hydration conformer correction",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
        (
            "confsolv_water_gsolv_std_teacher",
            "dispersion of conformer solvation energies",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
        (
            "confsolv_water_response_mean_teacher",
            "mean conformer solvent response",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
        (
            "confsolv_water_response_std_teacher",
            "dispersion of conformer solvent response",
            "ConfSolv H2O",
            "kcal mol-1",
            "LightGBM",
            39878,
            "direct surrogate",
        ),
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "prior_name",
            "physical_meaning",
            "source",
            "target_units",
            "model_class",
            "training_size",
            "transformation",
        ],
    )
    frame.insert(0, "prior", range(1, 16))
    frame["inference_computation"] = "structure-derived surrogate prediction; no simulation"
    return frame


def experiment_ledger(metrics: dict) -> pd.DataFrame:
    m = metrics["methods"]
    a = metrics["alternative_supervision"]
    ml = metrics["multilambda"]["method_mae_kcal_mol"]
    rows = [
        (
            "Previous structure-only",
            "Can molecular structure alone recover hydration free energy?",
            "RDKit descriptors + Morgan",
            "experimental hydration",
            "ExtraTrees",
            "structure only",
            "five-fold random OOF",
            "molecule-disjoint outer folds",
            m["previous_structure_only"]["mae_kcal_mol"],
            "baseline",
            "superseded",
            "results/predictions/headline_oof.parquet",
        ),
        (
            "Narrow response",
            "Do compact solvation-response coordinates improve the endpoint?",
            "CombiSolv-QM, Abraham, OpenFF, GBn2",
            "response priors + experimental hydration",
            "surrogate stack + ExtraTrees",
            "structure only",
            "five-fold random OOF",
            "external benchmark connectivity removed",
            m["narrow_response"]["mae_kcal_mol"],
            "improved",
            "final-model component",
            "results/predictions/headline_oof.parquet",
        ),
        (
            "SMD-water response",
            "Does aligned water response close the remaining gap?",
            "MolSolv SMD(water)",
            "SMD response + experimental hydration",
            "D-MPNN surrogate + ExtraTrees",
            "structure only",
            "five-fold random OOF",
            "85 benchmark connectivities removed",
            m["smd_water"]["mae_kcal_mol"],
            "improved",
            "final-model component",
            "results/predictions/headline_oof.parquet",
        ),
        (
            "ConfSolv response hierarchy",
            "Does conformational solvent response add transferable information?",
            "ConfSolv H2O",
            "six response summaries + experimental hydration",
            "LightGBM surrogates + ExtraTrees",
            "structure only",
            "five-fold random OOF",
            "13 benchmark connectivities removed",
            m["smd_confsolv_fixed"]["mae_kcal_mol"],
            "best fixed point estimate",
            "final-model component",
            "results/predictions/headline_oof.parquet",
        ),
        (
            "Nested feature-block selection",
            "Does feature-block selection preserve the result when nested?",
            "all retained response blocks",
            "experimental hydration",
            "nested ExtraTrees",
            "structure only",
            "nested five-fold OOF",
            "selection restricted to outer training data",
            m["nested_selection"]["mae_kcal_mol"],
            "PIMD8-level point estimate",
            "evaluation only",
            "results/predictions/headline_oof.parquet",
        ),
        (
            "Family holdout",
            "How does SolvAI transfer to an unseen functional family?",
            "final SolvAI",
            "experimental hydration",
            "frozen stack",
            "structure only",
            "GroupKFold family holdout",
            "families confined to one outer fold",
            m["family_holdout"]["mae_kcal_mol"],
            "harder chemical extrapolation",
            "evaluation only",
            "results/predictions/hard_holdout_oof.parquet",
        ),
        (
            "Scaffold holdout",
            "How does SolvAI transfer to unseen scaffolds?",
            "final SolvAI",
            "experimental hydration",
            "frozen stack",
            "structure only",
            "scaffold holdout",
            "scaffolds confined to one outer fold",
            m["scaffold_holdout"]["mae_kcal_mol"],
            "harder chemical extrapolation",
            "evaluation only",
            "results/predictions/hard_holdout_oof.parquet",
        ),
    ]
    for key, value in a.items():
        rows.append(
            (
                key.replace("_", " "),
                "Can this alternative physical representation improve transfer?",
                key,
                "experimental hydration",
                "documented frozen screen",
                "structure only",
                "matched OOF screen",
                "benchmark labels absent from held-out fit",
                float(value),
                "non-improving screen",
                "not retained",
                f"results/ablations/{key}.csv",
            )
        )
    for name, value in ml.items():
        rows.append(
            (
                name,
                "Can structure-derived lambda response improve the endpoint?",
                "short PIMD2 lambda response",
                "response hierarchy or integrated free energy",
                "masked multi-task student",
                "structure only",
                "matched random OOF",
                "held-out response labels absent",
                float(value),
                "response-surrogate bottleneck",
                "not retained",
                "results/ablations/multilambda_metrics.csv",
            )
        )
    return pd.DataFrame(
        rows,
        columns=[
            "experiment",
            "hypothesis",
            "physical_source",
            "target",
            "model",
            "inference_requirement",
            "evaluation_regime",
            "leakage_rule",
            "mae_kcal_mol",
            "result_interpretation",
            "role_in_final_model",
            "artifact_path",
        ],
    )


def main() -> None:
    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    methods, repeats = metrics["methods"], metrics["repeated_splits"]
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

    for folder in (TABLES, SI_TABLES, SUPP_DATA, ED):
        folder.mkdir(parents=True, exist_ok=True)
    macros = {
        "BenchmarkN": metrics["benchmark"]["molecules"],
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
    (TABLES / "metrics_macros.tex").write_text(
        "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items())
    )

    comparison = pd.DataFrame(
        [
            ("Classical ARROW", "yes", methods["classical_arrow"]["mae_kcal_mol"], "—", "—"),
            ("ARROW/PIMD8", "yes", methods["arrow_pimd8"]["mae_kcal_mol"], "—", "—"),
            (
                "Previous structure-only",
                "no",
                methods["previous_structure_only"]["mae_kcal_mol"],
                "—",
                "—",
            ),
            ("Narrow response", "no", methods["narrow_response"]["mae_kcal_mol"], "—", "—"),
            ("+ SMD(water)", "no", methods["smd_water"]["mae_kcal_mol"], "—", "—"),
            (
                "SolvAI",
                "no",
                methods["smd_confsolv_fixed"]["mae_kcal_mol"],
                methods["family_holdout"]["mae_kcal_mol"],
                methods["scaffold_holdout"]["mae_kcal_mol"],
            ),
            ("Nested selection", "no", methods["nested_selection"]["mae_kcal_mol"], "—", "—"),
        ],
        columns=[
            "Method",
            "Simulation at inference",
            "Random OOF MAE",
            "Family MAE",
            "Scaffold MAE",
        ],
    )
    comparison.to_csv(TABLES / "main_comparison.csv", index=False)
    comparison_fmt = comparison.copy()
    for c in ["Random OOF MAE", "Family MAE", "Scaffold MAE"]:
        comparison_fmt[c] = comparison_fmt[c].map(
            lambda x: f"{x:.3f}" if isinstance(x, float) else x
        )
    latex_table(comparison_fmt, TABLES / "main_comparison.tex", "lcccc")
    comparison.to_csv(ED / "ED_Table1.csv", index=False)
    latex_table(comparison_fmt, ED / "ED_Table1.tex", "lcccc")

    priors = response_priors()
    priors.to_csv(SUPP_DATA / "Supplementary_Data_4_teacher_priors.csv", index=False)
    latex_table(
        priors[
            ["prior", "prior_name", "physical_meaning", "source", "model_class", "training_size"]
        ],
        SI_TABLES / "response_priors.tex",
        "rllllr",
    )

    sources = pd.read_csv(ROOT / "data/manifests/training_source_manifest.csv")
    sources.to_csv(SUPP_DATA / "Supplementary_Data_4_teacher_sources.csv", index=False)
    latex_table(
        sources[
            [
                "source",
                "scientific_role",
                "original_records",
                "filtered_records",
                "benchmark_overlaps_removed",
                "license",
            ]
        ].fillna("—"),
        SI_TABLES / "source_provenance.tex",
        "llllrr",
    )

    endpoint = pd.DataFrame(
        [
            ("FreeSolv only", 69, "legacy benchmark-disjoint merge", "experimental hydration"),
            (
                "CombiSolv-Exp only",
                588,
                "legacy benchmark-disjoint merge",
                "experimental solvation in water",
            ),
            (
                "FreeSolv + CombiSolv-Exp identity",
                490,
                "one connectivity-level merged record",
                "duplicate-resolved experimental hydration",
            ),
            (
                "SoluteML dGsolvDB3",
                133,
                "at least two source measurements",
                "consensus experimental hydration",
            ),
        ],
        columns=["source_group", "selected_records", "selection_rule", "endpoint_target"],
    )
    endpoint.to_csv(SUPP_DATA / "Supplementary_Data_4_endpoint_sources.csv", index=False)
    latex_table(endpoint, SI_TABLES / "endpoint_sources.tex", "lrlp{4.2cm}")

    repeat_frame = pd.read_csv(ROOT / "results/robustness/repeated_metrics.csv")
    repeat_selected = repeat_frame.loc[
        repeat_frame.method.isin(["Nested selection", "narrow response + SMD + ConfSolv response"])
    ].copy()
    repeat_selected["mae"] = repeat_selected.mae.map(lambda x: f"{x:.5f}")
    latex_table(repeat_selected, SI_TABLES / "repeat_values.tex", "rrlrr")
    family = pd.DataFrame(metrics["chemistry_family"])
    family["mae_kcal_mol"] = family.mae_kcal_mol.map(lambda x: f"{x:.3f}")
    latex_table(family, SI_TABLES / "family_errors.tex", "lrr")

    manifest = json.loads((ROOT / "models/final/manifest.json").read_text())["model_files"]
    artifacts = pd.DataFrame([{"file": k, **v} for k, v in manifest.items()])
    latex_table(artifacts, SI_TABLES / "artifact_manifest.tex", "lrl")
    runtime = json.loads((ROOT / "results/runtime/runtime_benchmark.json").read_text())
    runtime_table = pd.DataFrame(
        [
            ("cold single molecule", f"{runtime['single_molecule']['cold_seconds']:.3f} s"),
            (
                "warm single molecule median",
                f"{runtime['single_molecule']['warm_median_seconds']:.3f} s",
            ),
            ("warm single molecule p95", f"{runtime['single_molecule']['warm_p95_seconds']:.3f} s"),
            ("batch 32 median", f"{runtime['batch']['median_seconds']:.3f} s"),
            ("batch 32 per molecule", f"{runtime['batch']['median_seconds_per_molecule']:.3f} s"),
            (
                "peak resident memory",
                f"{runtime['single_molecule']['peak_rss_kib'] / 1024:.1f} MiB",
            ),
        ],
        columns=["measurement", "value"],
    )
    latex_table(runtime_table, SI_TABLES / "runtime.tex", "ll")

    ledger = experiment_ledger(metrics)
    ledger.to_csv(SUPP_DATA / "Supplementary_Data_1_experiment_ledger.csv", index=False)
    workbook(SUPP_DATA / "Supplementary_Data_1_experiment_ledger.xlsx", {"experiments": ledger})

    headline = pd.read_parquet(ROOT / "results/predictions/headline_oof.parquet")
    base = headline[
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            "fold",
            "y_true",
            "delta_g_pimd8",
            "delta_g_classical_arrow",
        ]
    ].drop_duplicates("molecule_id")
    wide = headline.pivot(index="molecule_id", columns="method", values="y_pred").reset_index()
    molecule_predictions = base.merge(wide, on="molecule_id", validate="one_to_one")
    molecule_predictions.to_csv(
        SUPP_DATA / "Supplementary_Data_2_molecule_predictions.csv", index=False
    )
    workbook(
        SUPP_DATA / "Supplementary_Data_2_molecule_predictions.xlsx",
        {"predictions": molecule_predictions},
    )

    benchmark = pd.read_parquet(ROOT / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[
        benchmark.solvent.eq("water"),
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            "fold_random",
            "fold_family",
            "fold_scaffold",
        ],
    ]
    repeated = pd.read_parquet(ROOT / "results/robustness/repeated_oof.parquet")
    repeated = repeated.loc[
        repeated.method.eq("narrow response + SMD + ConfSolv response"),
        ["molecule_id", "repeat", "split_seed", "outer_fold"],
    ]
    repeat_wide = (
        repeated.pivot(index="molecule_id", columns="repeat", values="outer_fold")
        .add_prefix("repeat_fold_")
        .reset_index()
    )
    assignments = benchmark.merge(repeat_wide, on="molecule_id", validate="one_to_one")
    assignments.to_csv(SUPP_DATA / "Supplementary_Data_3_split_assignments.csv", index=False)
    workbook(
        SUPP_DATA / "Supplementary_Data_3_split_assignments.xlsx", {"assignments": assignments}
    )

    workbook(
        SUPP_DATA / "Supplementary_Data_4_teacher_manifests.xlsx",
        {"response_priors": priors, "teacher_sources": sources, "endpoint_sources": endpoint},
    )
    print("Rendered manuscript tables and four machine-readable Supplementary Data packages.")


if __name__ == "__main__":
    main()
