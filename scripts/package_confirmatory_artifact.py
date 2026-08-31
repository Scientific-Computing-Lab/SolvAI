#!/usr/bin/env python3
"""Package the standardized-exclusion SolvAI artifact without model selection."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from confirmatory_common import MODEL_SEEDS, endpoint_model, load_confirmatory_data

RELEASE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = RELEASE_ROOT.parents[1]
MODEL_DIR = RELEASE_ROOT / "models" / "final"
CONFIRMATORY = RELEASE_ROOT / "results" / "confirmatory"


def main() -> None:
    overrides = {
        "combisolv_qm": CONFIRMATORY / "teacher_refits/combisolv_qm/teacher_predictions.parquet",
        "molsolv_smd": CONFIRMATORY / "teacher_refits/molsolv_smd/teacher_predictions.parquet",
        "confsolv": CONFIRMATORY / "teacher_refits/confsolv/teacher_predictions.parquet",
    }
    data = load_confirmatory_data(WORKSPACE_ROOT, overrides)
    benchmark_x, public_x = data.feature_sets["F_full_solvai"]
    x_fit = np.vstack([public_x, benchmark_x])
    y_fit = np.concatenate(
        [
            data.public.delta_g_exp.to_numpy(dtype=float),
            data.benchmark.delta_g_exp.to_numpy(dtype=float),
        ]
    )
    weights = np.concatenate([np.ones(len(data.public)), np.full(len(data.benchmark), 3.0)])
    models = []
    for seed in MODEL_SEEDS:
        model = endpoint_model(seed)
        model.fit(x_fit, y_fit, extratreesregressor__sample_weight=weights)
        models.append(model)

    descriptor_columns = [
        column
        for column in pd.read_parquet(
            WORKSPACE_ROOT / "data/processed/rdkit_morgan_features.parquet"
        )
        if column.startswith(("rdkit__", "morgan2__"))
    ]
    if len(descriptor_columns) != 2265 or benchmark_x.shape[1] != 2280:
        raise AssertionError("Confirmatory artifact feature schema changed")
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": descriptor_columns,
            "teacher_columns": data.response_names,
            "feature_block": "full retained SolvAI (standardized-exclusion teachers)",
            "seeds": list(MODEL_SEEDS),
            "trees": 360,
        },
        MODEL_DIR / "head.joblib",
        compress=3,
    )

    replacements = {
        "combisolv_qm.pt": CONFIRMATORY / "teacher_refits/combisolv_qm/model/model_0/best.pt",
        "molsolv_smd.pt": CONFIRMATORY / "teacher_refits/molsolv_smd/model/model_0/best.pt",
        "confsolv_teacher.joblib": CONFIRMATORY / "teacher_refits/confsolv/lightgbm_models.joblib",
    }
    for destination, source in replacements.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, MODEL_DIR / destination)

    summary = json.loads((CONFIRMATORY / "confirmatory_summary.json").read_text())
    external_metrics = pd.read_csv(
        RELEASE_ROOT / "results/tier_a_external/evaluation/tier_a_external_metrics.csv"
    )
    endpoint_external = external_metrics.loc[
        external_metrics.cohort.eq("endpoint_disjoint")
        & external_metrics.method.eq("full_solvai")
    ].iloc[0]
    strict_external = external_metrics.loc[
        external_metrics.cohort.eq("strict_response_source_disjoint")
        & external_metrics.method.eq("full_solvai")
    ].iloc[0]
    card = {
        "artifact": "models/final/head.joblib",
        "input": "SMILES only",
        "output": "hydration free energy in kcal/mol",
        "input_canonicalization": "RDKit canonical isomeric SMILES",
        "inference_simulation": False,
        "inference_forbidden_inputs": [
            "MD",
            "PIMD",
            "ARROW trajectory",
            "physics probe",
            "experimental label",
        ],
        "selected_feature_block": "15 predicted solvation-response priors",
        "confirmatory_fixed_oof_mae": summary["primary"]["full_solvai_mae"],
        "matched_structure_only_oof_mae": summary["primary"]["matched_structure_only_mae"],
        "five_repeat_mean_mae": summary["repeats"]["full_solvai_mean"],
        "five_repeat_sd_mae": summary["repeats"]["full_solvai_sd"],
        "zero_arrow_label_transfer_mae": summary["zero_arrow_transfer"]["F_full_solvai"],
        "tier_a_endpoint_disjoint_n": int(endpoint_external.n),
        "tier_a_endpoint_disjoint_mae": float(endpoint_external.mae),
        "tier_a_strict_source_disjoint_n": int(strict_external.n),
        "tier_a_strict_source_disjoint_mae": float(strict_external.mae),
        "global_family_mae": next(
            row["mae"]
            for row in summary["global_separation"]
            if row["regime"] == "global_family" and row["method"] == "F_full_solvai"
        ),
        "global_scaffold_mae": next(
            row["mae"]
            for row in summary["global_separation"]
            if row["regime"] == "global_scaffold" and row["method"] == "F_full_solvai"
        ),
        "robust_below_0_20": False,
        "scikit_learn_version": sklearn.__version__,
        "benchmark_labels_in_post_evaluation_refit": len(data.benchmark),
        "public_labels_in_refit": len(data.public),
        "pimd8_labels_in_selected_artifact": 0,
        "calibrated_applicability_score": False,
        "pimd8_note": (
            "PIMD-derived supervision was evaluated but not retained. ARROW/PIMD8 "
            "is an accuracy comparator, not an input to the released model."
        ),
        "standardized_teacher_exclusions": summary["standardized_teacher_exclusions"],
        "teacher_artifacts": [
            "combisolv_qm.pt",
            "molsolv_smd.pt",
            "abraham_teacher.joblib",
            "openff_teacher.joblib",
            "implicit_teacher.joblib",
            "confsolv_teacher.joblib",
        ],
        "evaluation_predictions": (
            "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
        ),
    }
    (MODEL_DIR / "model_card.json").write_text(json.dumps(card, indent=2) + "\n")
    (MODEL_DIR / "MODEL_CARD.md").write_text(
        f"""# SolvAI model card

SolvAI predicts neutral-molecule hydration free energy from a SMILES string. The
released stack computes deterministic molecular descriptors and 15 response priors
with structure surrogates, then applies an ensemble of three ExtraTrees endpoint
models. No MD, PIMD, ARROW trajectory, probe calculation or experimental lookup is
performed at inference.

RDKit parses each query and emits canonical isomeric SMILES before feature
generation. This is not a claim of invariance across tautomers, protonation states or
salt forms.

## Validated performance

- Confirmatory fixed five-fold OOF MAE: {card["confirmatory_fixed_oof_mae"]:.5f} kcal/mol
- Matched no-prior endpoint MAE: {card["matched_structure_only_oof_mae"]:.5f} kcal/mol
- Five-partition mean: {card["five_repeat_mean_mae"]:.5f} ± {card["five_repeat_sd_mae"]:.5f} kcal/mol
- Zero-ARROW-label transfer MAE: {card["zero_arrow_label_transfer_mae"]:.5f} kcal/mol
- Tier-A endpoint-disjoint MAE (N={card["tier_a_endpoint_disjoint_n"]}): {card["tier_a_endpoint_disjoint_mae"]:.5f} kcal/mol
- Tier-A strict source-disjoint MAE (N={card["tier_a_strict_source_disjoint_n"]}): {card["tier_a_strict_source_disjoint_mae"]:.5f} kcal/mol
- Global family / scaffold MAE: {card["global_family_mae"]:.5f} / {card["global_scaffold_mae"]:.5f} kcal/mol

The reference domain is neutral organic hydration chemistry. The model is not
validated for ions, salts, metals, proteins or broad chemical extrapolation. The
fixed point estimate is on the ARROW/PIMD8 accuracy scale, but robust sub-0.20
performance and superiority over PIMD8 are not claimed. The returned ensemble spread
is not a calibrated per-query applicability or reliability score.

## Training and inference boundary

Physical calculations supplied training targets for the response surrogates. The
deployed model receives structure only. PIMD-derived candidate features were tested
and not retained; PIMD8 is used solely as an accuracy comparator.
"""
    )


if __name__ == "__main__":
    main()
