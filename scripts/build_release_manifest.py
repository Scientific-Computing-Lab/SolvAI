#!/usr/bin/env python3
"""Hash the canonical scientific and publication artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = [
    "data/benchmark/arrow_solvation_master.parquet",
    "data/manifests/endpoint_label_manifest.json",
    "data/manifests/training_source_manifest.json",
    "results/paper_metrics.json",
    "results/paper_metrics.csv",
    "results/final_metrics.json",
    "results/model_comparison.csv",
    "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet",
    "results/confirmatory/standardized_exclusion_endpoint_shuffle_predictions.parquet",
    "results/confirmatory/standardized_exclusion_global_separation_predictions.parquet",
    "results/confirmatory/zero_arrow_transfer_predictions.parquet",
    "results/confirmatory/confirmatory_summary.json",
    "results/michael_30aug_sensitivity/weight1_predictions.parquet",
    "results/michael_30aug_sensitivity/weight1_metrics.csv",
    "results/michael_30aug_sensitivity/weight1_paired_comparison.csv",
    "results/michael_30aug_sensitivity/weight1_metadata.json",
    "results/tier_a_external/qualification/tier_a_qualification_audit.csv",
    "results/tier_a_external/qualification/tier_a_qualification_audit.parquet",
    "results/tier_a_external/qualification/tier_a_qualification_summary.json",
    "results/tier_a_external/qualification/tier_a_identity_matches.csv",
    "results/tier_a_external/qualification/tier_a_exclusions.csv",
    "results/tier_a_external/qualification/tier_a_endpoint_disjoint.csv",
    "results/tier_a_external/qualification/tier_a_endpoint_disjoint.parquet",
    "results/tier_a_external/qualification/tier_a_strict_response_source_disjoint.csv",
    "results/tier_a_external/qualification/tier_a_strict_response_source_disjoint.parquet",
    "results/tier_a_external/qualification/tier_a_teacher_exposure_summary.csv",
    "results/tier_a_external/evaluation/tier_a_external_predictions.csv",
    "results/tier_a_external/evaluation/tier_a_external_predictions.parquet",
    "results/tier_a_external/evaluation/tier_a_external_response_features.parquet",
    "results/tier_a_external/evaluation/tier_a_external_metrics.csv",
    "results/tier_a_external/evaluation/tier_a_external_paired_comparisons.csv",
    "results/tier_a_external/evaluation/tier_a_external_descriptive_strata.csv",
    "results/tier_a_external/evaluation/tier_a_external_metadata.json",
    "audits/leakage_audit.json",
    "audits/confirmatory/chemical_distance_audit.json",
    "audits/confirmatory/standardized_exclusion_refit_verification.json",
    "audits/artifact_audit.json",
    "audits/claim_red_team.json",
    "audits/security_audit.json",
    "models/final/manifest.json",
    "results/runtime/runtime_benchmark.json",
    "release/CONFIRMATORY_FREEZE.md",
    "release/MICHAEL_30AUG_SENSITIVITY_FREEZE.md",
    "release/TIER_A_EXTERNAL_VALIDATION_FREEZE.md",
    "reports/CONFIRMATORY_ANALYSIS.md",
    "reports/MICHAEL_30AUG_WEIGHT_SENSITIVITY.md",
    "reports/TIER_A_EXTERNAL_VALIDATION.md",
    "reports/PAPER_FREEZE.md",
    "README.md",
    "paper/main.tex",
    "paper/references.bib",
    "paper/main.pdf",
    "paper/supplementary/supplementary.pdf",
    "paper/supplementary/supplementary.tex",
    "paper/supplementary_data/Supplementary_Data_1_experiment_ledger.xlsx",
    "paper/supplementary_data/Supplementary_Data_2_molecule_predictions.xlsx",
    "paper/supplementary_data/Supplementary_Data_3_split_assignments.xlsx",
    "paper/supplementary_data/Supplementary_Data_4_teacher_manifests.xlsx",
    "paper/supplementary_data/Supplementary_Data_5_tier_a_external_validation.xlsx",
    "paper/figures/main/fig1_concept.pdf",
    "paper/figures/main/fig2_headline.pdf",
    "paper/figures/main/fig3_transfer.pdf",
    "paper/figures/main/fig4_frontier.pdf",
    "paper/extended_data/ED_Fig1_residuals.pdf",
    "paper/extended_data/ED_Fig2_provenance.pdf",
    "paper/extended_data/ED_Fig3_alternatives.pdf",
    "paper/extended_data/ED_Fig4_selective_pimd.pdf",
    "paper/extended_data/ED_Fig5_lambda_response.pdf",
    "paper/extended_data/ED_Fig6_extrapolation.pdf",
    "paper/extended_data/ED_Table1.pdf",
    "paper/review_combined.pdf",
    "paper/supplementary_data/Supplementary_Data_1_experiment_ledger.csv",
    "paper/supplementary_data/Supplementary_Data_2_molecule_predictions.csv",
    "paper/supplementary_data/Supplementary_Data_3_split_assignments.csv",
    "paper/supplementary_data/Supplementary_Data_4_teacher_priors.csv",
    "paper/supplementary_data/Supplementary_Data_5_tier_a_predictions.csv",
    "paper/supplementary_data/Supplementary_Data_5_tier_a_qualification.csv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = {}
    for relative in TARGETS:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        files[relative] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {"schema_version": 1, "files": files}
    output = ROOT / "release_manifest.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Hashed {len(files)} canonical artifacts")


if __name__ == "__main__":
    main()
