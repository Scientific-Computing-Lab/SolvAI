"""Fit and package the frozen SMD+ConfSolv structure-only student.

The reported score for this model comes from strict OOF predictions.  This
script is only the post-evaluation refit on all available labels for deployment.
Every packaged teacher maps molecular structure to a training-derived physical
quantity; no trajectory or measured property is requested at inference.
"""

from __future__ import annotations

import argparse
import json
import shutil

import joblib
import pandas as pd
import sklearn
from arrow_distill.data import ROOT
from run_nested_smd_teacher_confirmation import fit_model, load_problem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", type=int, default=360)
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 29, 47])
    args = parser.parse_args()

    benchmark, public, feature_sets, source_mask = load_problem(
        no_static_arrow=True,
        confsolv_response=True,
    )
    feature_name = "narrow response + SMD + ConfSolv response"
    benchmark_x, public_x = feature_sets[feature_name]
    benchmark_y = benchmark.delta_g_exp.to_numpy(dtype=float)
    public_y = public.delta_g_exp.to_numpy(dtype=float)[source_mask]
    models = [
        fit_model(
            public_x[source_mask],
            public_y,
            benchmark_x,
            benchmark_y,
            seed,
            args.trees,
        )
        for seed in args.seeds
    ]

    base = pd.read_parquet(ROOT / "data/processed/rdkit_morgan_features.parquet")
    descriptor_columns = [column for column in base if column.startswith(("rdkit__", "morgan2__"))]
    teacher_columns = [
        "combisolv_qm_teacher",
        "abraham_e_teacher",
        "abraham_s_teacher",
        "abraham_a_teacher",
        "abraham_b_teacher",
        "abraham_l_teacher",
        "openff_corrected_teacher",
        "gbn2_corrected_teacher",
        "molsolv_smd_teacher",
        "confsolv_gas_conformer_correction_teacher",
        "confsolv_solution_conformer_correction_teacher",
        "confsolv_hydration_conformer_correction_teacher",
        "confsolv_water_gsolv_std_teacher",
        "confsolv_water_response_mean_teacher",
        "confsolv_water_response_std_teacher",
    ]
    if benchmark_x.shape[1] != len(descriptor_columns) + len(teacher_columns):
        raise AssertionError(
            f"Feature schema mismatch: matrix={benchmark_x.shape[1]}, "
            f"schema={len(descriptor_columns) + len(teacher_columns)}"
        )

    output_dir = ROOT / "models/final_structure_only"
    output_dir.mkdir(parents=True, exist_ok=True)
    head_path = output_dir / "head.joblib"
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": descriptor_columns,
            "teacher_columns": teacher_columns,
            "feature_block": feature_name,
            "seeds": args.seeds,
            "trees": args.trees,
        },
        head_path,
        compress=3,
    )

    dependencies = {
        "combisolv_qm.pt": ROOT / "results/combisolv_qm_pretraining/model/model_0/best.pt",
        "molsolv_smd.pt": ROOT / "results/molsolv_smd_pretraining/model/model_0/best.pt",
        "abraham_teacher.joblib": ROOT / "models/soluteml_abraham_teacher/models.joblib",
        "openff_teacher.joblib": ROOT / "models/openff_alchemical_teacher/models.joblib",
        "implicit_teacher.joblib": ROOT / "models/implicit_solvent_teacher/models.joblib",
        "confsolv_teacher.joblib": ROOT / "models/confsolv_water_teacher/lightgbm_models.joblib",
    }
    missing = [str(path) for path in dependencies.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Build the serialized teacher artifacts before packaging: " + ", ".join(missing)
        )
    for destination, source in dependencies.items():
        shutil.copy2(source, output_dir / destination)

    confirmation = pd.read_csv(ROOT / "results/nested_smd_confsolv_structure_confirmation.csv")
    hard = pd.read_csv(ROOT / "results/smd_confsolv_response_structure_hard.csv")
    repeated = pd.read_csv(
        ROOT / "results/repeated_nested_smd_confsolv_confirmation_summary.csv"
    ).set_index("method")
    metadata = {
        "artifact": str(head_path.relative_to(ROOT)),
        "input": "SMILES only",
        "output": "hydration free energy in kcal/mol",
        "inference_simulation": False,
        "inference_forbidden_inputs": [
            "MD",
            "PIMD",
            "ARROW trajectory",
            "physics probe",
            "experimental label",
        ],
        "selected_feature_block": feature_name,
        "nested_random_oof_mae": float(
            confirmation.loc[
                confirmation.method.eq("Nested narrow SMD teacher selection"), "mae"
            ].iloc[0]
        ),
        "fixed_random_oof_mae": float(
            confirmation.loc[
                confirmation.method.eq("Fixed narrow response + SMD + ConfSolv response"),
                "mae",
            ].iloc[0]
        ),
        "fixed_family_holdout_mae": float(
            hard.loc[hard.regime.eq("family_holdout"), "mae"].iloc[0]
        ),
        "fixed_scaffold_holdout_mae": float(
            hard.loc[hard.regime.eq("scaffold_holdout"), "mae"].iloc[0]
        ),
        "fixed_repeat_mean_mae": float(
            repeated.loc["narrow response + SMD + ConfSolv response", "mean_mae"]
        ),
        "nested_repeat_mean_mae": float(repeated.loc["Nested selection", "mean_mae"]),
        "robust_below_0_20": False,
        "scikit_learn_version": sklearn.__version__,
        "benchmark_labels_in_post_evaluation_refit": len(benchmark),
        "public_labels_in_refit": int(source_mask.sum()),
        "pimd8_labels_in_selected_artifact": 0,
        "pimd8_note": (
            "85 benchmark PIMD8 labels were evaluated as training-only candidate "
            "supervision, but that block worsened OOF validation and was not selected."
        ),
        "teacher_artifacts": list(dependencies),
        "evaluation_predictions": "results/nested_smd_confsolv_structure_confirmation.parquet",
    }
    (output_dir / "model_card.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
