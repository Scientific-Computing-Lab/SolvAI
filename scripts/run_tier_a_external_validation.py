#!/usr/bin/env python3
"""Run the prospectively frozen Tier-A matched external validation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from confirmatory_common import (
    BOOTSTRAP_SEED,
    MODEL_SEEDS,
    bootstrap_difference,
    endpoint_model,
    load_confirmatory_data,
    metric_record,
)
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

from solv_ai.features import canonicalize, descriptor_frame
from solv_ai.teachers import response_features

RELEASE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = RELEASE_ROOT.parents[1]
QUALIFICATION = RELEASE_ROOT / "results/tier_a_external/qualification"
OUT = RELEASE_ROOT / "results/tier_a_external/evaluation"
MODEL_DIR = RELEASE_ROOT / "models/final"
PROTOCOL_COMMIT = "b01351a"
BOOTSTRAP_DRAWS = 100_000

EXPECTED_INPUT_HASHES = {
    QUALIFICATION
    / "tier_a_endpoint_disjoint.csv": "967a9794a5d0e3f131dc5bb6921fecd234809ca21143a5a3ce9ed7895d18b273",
    QUALIFICATION
    / "tier_a_strict_response_source_disjoint.csv": "486040f54082f606ca21f461eb76075648b885c712823c72f7f69317960463bf",
    WORKSPACE_ROOT
    / "data/processed/arrow_solvation_master.parquet": "2b7928f162d094e7ee10d197e66636ba4ae09b0f76d626136b79c0975d3b0310",
    WORKSPACE_ROOT
    / "data/processed/expanded_public_hydration_nonbenchmark.parquet": "603ed02b6be25d9a3057e321f2c6ea135b012666cfdb8a1b160e37f347951ec4",
    MODEL_DIR
    / "combisolv_qm.pt": "2950be1232e2f868136e3388cf139e671d4196922ece0d3c7043d9dadece8b40",
    MODEL_DIR
    / "molsolv_smd.pt": "35601030f030af72a9e6990bf8003df9d481858b27994b3fc70c81d605d130be",
    MODEL_DIR
    / "abraham_teacher.joblib": "8c8a6365e7162dbe6ae445ce4004928ed86b2669d3e7a70d3d054f143ca97554",
    MODEL_DIR
    / "openff_teacher.joblib": "55b7d31a293dd61fe57f5505e015775d7d86f53dca69980d398843706d58bb9b",
    MODEL_DIR
    / "implicit_teacher.joblib": "484d5e0c1efa759c8e7acbfcbbd87e218b93ae239e7f35303fe4cd7e310910f8",
    MODEL_DIR
    / "confsolv_teacher.joblib": "a71487e53d5085e897dad02e62e95562933bc47262145fbf44286aae04b70b0d",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    """Describe an input without embedding a release-host home directory."""
    try:
        return f"repository/{path.relative_to(RELEASE_ROOT)}"
    except ValueError:
        pass
    try:
        return f"workspace/{path.relative_to(WORKSPACE_ROOT)}"
    except ValueError:
        return path.name


def teacher_overrides() -> dict[str, Path]:
    root = RELEASE_ROOT / "results/confirmatory/teacher_refits"
    return {
        "combisolv_qm": root / "combisolv_qm/teacher_predictions.parquet",
        "molsolv_smd": root / "molsolv_smd/teacher_predictions.parquet",
        "confsolv": root / "confsolv/teacher_predictions.parquet",
    }


def predict_ensemble(
    train_x: np.ndarray,
    train_y: np.ndarray,
    weights: np.ndarray,
    test_x: np.ndarray,
) -> tuple[np.ndarray, list[object]]:
    predictions = []
    models = []
    for seed in MODEL_SEEDS:
        model = endpoint_model(seed)
        model.fit(train_x, train_y, extratreesregressor__sample_weight=weights)
        models.append(model)
        predictions.append(model.predict(test_x))
    return np.mean(predictions, axis=0), models


def nearest_endpoint_similarity(
    training_smiles: list[str], candidate_smiles: list[str]
) -> tuple[np.ndarray, list[str]]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    training = [
        generator.GetFingerprint(Chem.MolFromSmiles(canonicalize(value)[0]))
        for value in training_smiles
    ]
    maximum = np.empty(len(candidate_smiles), dtype=float)
    nearest: list[str] = []
    for index, value in enumerate(candidate_smiles):
        fingerprint = generator.GetFingerprint(Chem.MolFromSmiles(canonicalize(value)[0]))
        similarities = np.asarray(
            DataStructs.BulkTanimotoSimilarity(fingerprint, training), dtype=float
        )
        nearest_index = int(np.argmax(similarities))
        maximum[index] = similarities[nearest_index]
        nearest.append(training_smiles[nearest_index])
    return maximum, nearest


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path, expected in EXPECTED_INPUT_HASHES.items():
        observed = sha256(path)
        if observed != expected:
            raise AssertionError(f"Frozen input changed for {path}: {observed}")

    external = pd.read_csv(QUALIFICATION / "tier_a_endpoint_disjoint.csv")
    strict_ids = set(
        pd.read_csv(QUALIFICATION / "tier_a_strict_response_source_disjoint.csv").candidate_id
    )
    if len(external) != 220 or len(strict_ids) != 97:
        raise AssertionError("Frozen Tier-A cohort counts changed")
    external["strict_response_source_disjoint"] = external.candidate_id.isin(strict_ids)
    canonical = [canonicalize(value)[0] for value in external.smiles.astype(str)]
    descriptors = descriptor_frame(canonical)
    head = joblib.load(MODEL_DIR / "head.joblib")
    descriptor_columns = list(head["descriptor_columns"])
    missing = set(descriptor_columns) - set(descriptors.columns)
    if missing:
        raise AssertionError(f"External descriptors are missing {len(missing)} columns")
    external_structure = descriptors[descriptor_columns].to_numpy(dtype=np.float32)
    external_responses = response_features(canonical, descriptors, MODEL_DIR).astype(np.float32)
    if external_structure.shape != (220, 2265) or external_responses.shape != (220, 15):
        raise AssertionError("Unexpected Tier-A feature shape")
    if not np.isfinite(external_responses).all():
        raise AssertionError("Non-finite Tier-A response prediction")

    data = load_confirmatory_data(WORKSPACE_ROOT, teacher_overrides())
    structure_benchmark, structure_public = data.feature_sets["A_structure_only"]
    full_benchmark, full_public = data.feature_sets["F_full_solvai"]
    train_y = np.concatenate(
        [
            data.public.delta_g_exp.to_numpy(dtype=float),
            data.benchmark.delta_g_exp.to_numpy(dtype=float),
        ]
    )
    weights = np.concatenate(
        [np.ones(len(data.public)), np.full(len(data.benchmark), 3.0)]
    )
    structure_train = np.vstack([structure_public, structure_benchmark])
    full_train = np.vstack([full_public, full_benchmark])
    structure_prediction, _ = predict_ensemble(
        structure_train, train_y, weights, external_structure
    )
    full_test = np.column_stack([external_structure, external_responses])
    full_prediction, _ = predict_ensemble(full_train, train_y, weights, full_test)

    released_prediction = np.mean(
        [model.predict(full_test) for model in head["models"]], axis=0
    )
    artifact_maximum_difference = float(np.max(np.abs(released_prediction - full_prediction)))
    if artifact_maximum_difference > 1e-10:
        raise AssertionError(
            "Refitted full endpoint does not reproduce the released endpoint: "
            f"max difference {artifact_maximum_difference}"
        )

    training_smiles = list(data.public.canonical_smiles.astype(str)) + list(
        data.benchmark.canonical_smiles.astype(str)
    )
    maximum_similarity, nearest_smiles = nearest_endpoint_similarity(training_smiles, canonical)
    predictions = external.copy()
    predictions["canonical_smiles_evaluated"] = canonical
    predictions["y_true"] = predictions.experimental_dg_kcal_mol.astype(float)
    predictions["structure_only_prediction"] = structure_prediction
    predictions["solvai_prediction"] = full_prediction
    predictions["structure_only_residual"] = predictions.y_true - structure_prediction
    predictions["solvai_residual"] = predictions.y_true - full_prediction
    predictions["structure_only_absolute_error"] = predictions.structure_only_residual.abs()
    predictions["solvai_absolute_error"] = predictions.solvai_residual.abs()
    predictions["solvai_improved"] = (
        predictions.solvai_absolute_error < predictions.structure_only_absolute_error
    )
    predictions["maximum_endpoint_training_morgan_tanimoto"] = maximum_similarity
    predictions["nearest_endpoint_training_smiles"] = nearest_smiles
    predictions.to_csv(OUT / "tier_a_external_predictions.csv", index=False)
    predictions.to_parquet(OUT / "tier_a_external_predictions.parquet", index=False)

    response_table = pd.DataFrame(external_responses, columns=data.response_names)
    response_table.insert(0, "canonical_smiles", canonical)
    response_table.insert(0, "candidate_id", external.candidate_id)
    response_table.to_parquet(OUT / "tier_a_external_response_features.parquet", index=False)

    metrics: list[dict[str, object]] = []
    comparisons: list[dict[str, object]] = []
    for cohort_name, mask in (
        ("endpoint_disjoint", np.ones(len(predictions), dtype=bool)),
        (
            "strict_response_source_disjoint",
            predictions.strict_response_source_disjoint.to_numpy(dtype=bool),
        ),
    ):
        subset = predictions.loc[mask]
        truth = subset.y_true.to_numpy(dtype=float)
        structure = subset.structure_only_prediction.to_numpy(dtype=float)
        solvai = subset.solvai_prediction.to_numpy(dtype=float)
        for method, values in (
            ("matched_structure_only", structure),
            ("full_solvai", solvai),
        ):
            metrics.append(
                {
                    "cohort": cohort_name,
                    "method": method,
                    "n": len(subset),
                    **metric_record(truth, values),
                }
            )
        comparisons.append(
            {
                "cohort": cohort_name,
                "n": len(subset),
                **bootstrap_difference(
                    truth,
                    solvai,
                    structure,
                    seed=BOOTSTRAP_SEED,
                    draws=BOOTSTRAP_DRAWS,
                ),
            }
        )
    metrics_frame = pd.DataFrame(metrics)
    comparisons_frame = pd.DataFrame(comparisons)
    metrics_frame.to_csv(OUT / "tier_a_external_metrics.csv", index=False)
    comparisons_frame.to_csv(OUT / "tier_a_external_paired_comparisons.csv", index=False)

    size_bins = pd.qcut(predictions.heavy_atom_count, q=4, duplicates="drop")
    similarity_bins = pd.qcut(
        predictions.maximum_endpoint_training_morgan_tanimoto, q=4, duplicates="drop"
    )
    descriptive_rows: list[dict[str, object]] = []
    for dimension, bins in (
        ("heavy_atom_count_quartile", size_bins),
        ("nearest_similarity_quartile", similarity_bins),
    ):
        for group, subset in predictions.groupby(bins, observed=True):
            descriptive_rows.append(
                {
                    "dimension": dimension,
                    "group": str(group),
                    "n": len(subset),
                    "structure_only_mae": float(subset.structure_only_absolute_error.mean()),
                    "solvai_mae": float(subset.solvai_absolute_error.mean()),
                    "paired_mae_difference": float(
                        (subset.solvai_absolute_error - subset.structure_only_absolute_error).mean()
                    ),
                }
            )
    pd.DataFrame(descriptive_rows).to_csv(
        OUT / "tier_a_external_descriptive_strata.csv", index=False
    )

    output_files = [
        "tier_a_external_predictions.csv",
        "tier_a_external_predictions.parquet",
        "tier_a_external_response_features.parquet",
        "tier_a_external_metrics.csv",
        "tier_a_external_paired_comparisons.csv",
        "tier_a_external_descriptive_strata.csv",
    ]
    metadata = {
        "protocol": "release/TIER_A_EXTERNAL_VALIDATION_FREEZE.md",
        "protocol_commit": PROTOCOL_COMMIT,
        "endpoint_training_rows": {"external": len(data.public), "arrow": len(data.benchmark)},
        "endpoint_weights": {"external": 1, "arrow": 3},
        "cohort_rows": {"endpoint_disjoint": 220, "strict_response_source_disjoint": 97},
        "model_seeds": list(MODEL_SEEDS),
        "bootstrap": {"draws": BOOTSTRAP_DRAWS, "seed": BOOTSTRAP_SEED},
        "released_artifact_maximum_prediction_difference": artifact_maximum_difference,
        "rdkit_version": rdBase.rdkitVersion,
        "input_sha256": {
            portable_path(path): value for path, value in EXPECTED_INPUT_HASHES.items()
        },
        "output_sha256": {name: sha256(OUT / name) for name in output_files},
        "metrics": metrics,
        "paired_comparisons": comparisons,
    }
    (OUT / "tier_a_external_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(metrics_frame.to_string(index=False))
    print(comparisons_frame.to_string(index=False))
    print(f"released artifact max prediction difference: {artifact_maximum_difference:.3e}")


if __name__ == "__main__":
    main()
