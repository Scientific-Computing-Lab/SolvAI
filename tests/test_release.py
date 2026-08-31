from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pytest

from solv_ai import predict_smiles
from solv_ai.cli import parse_smiles
from solv_ai.features import canonicalize, descriptor_frame
from solv_ai.paper_metrics import compute_paper_metrics

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def smoke_prediction() -> tuple[np.ndarray, np.ndarray]:
    return predict_smiles(["CCO", "c1ccccc1"])


def test_canonicalization_is_deterministic() -> None:
    canonical_a, key_a, connectivity_a = canonicalize("OCC")
    canonical_b, key_b, connectivity_b = canonicalize("CCO")
    assert canonical_a == canonical_b == "CCO"
    assert key_a == key_b
    assert connectivity_a == connectivity_b == key_a.split("-")[0]


def test_invalid_smiles_is_rejected() -> None:
    with pytest.raises(ValueError, match="could not parse"):
        canonicalize("not-a-smiles")


def test_empty_inference_is_rejected() -> None:
    with pytest.raises(ValueError, match="At least one"):
        predict_smiles([])


def test_cli_accepts_documented_predict_subcommand() -> None:
    assert parse_smiles(["predict", "CCO"]) == ["CCO"]
    assert parse_smiles(["CCO"]) == ["CCO"]


def test_descriptor_schema_matches_artifact() -> None:
    artifact = joblib.load(ROOT / "models/final/head.joblib")
    frame = descriptor_frame(["CCO"])
    assert len(artifact["descriptor_columns"]) == 2265
    assert set(artifact["descriptor_columns"]).issubset(frame.columns)


def test_artifact_contains_no_prohibited_feature() -> None:
    artifact = joblib.load(ROOT / "models/final/head.joblib")
    names = [*artifact["descriptor_columns"], *artifact["teacher_columns"]]
    prohibited = ("delta_g_exp", "pimd", "trajectory", "probe", "fold_", "scaffold")
    assert not [name for name in names if any(term in name.lower() for term in prohibited)]


def test_smiles_only_smoke_predictions(smoke_prediction: tuple[np.ndarray, np.ndarray]) -> None:
    prediction, spread = smoke_prediction
    manifest = json.loads((ROOT / "models/final/manifest.json").read_text())
    expected_prediction = [row["prediction_kcal_mol"] for row in manifest["smoke_test"]]
    expected_spread = [row["ensemble_spread_kcal_mol"] for row in manifest["smoke_test"]]
    np.testing.assert_allclose(prediction, expected_prediction, atol=2e-6, rtol=0)
    np.testing.assert_allclose(spread, expected_spread, atol=2e-6, rtol=0)


def test_paper_metrics_recompute_without_drift() -> None:
    metrics, table = compute_paper_metrics(ROOT)
    assert metrics["benchmark"]["molecules"] == 85
    assert metrics["methods"]["full_solvai"]["mae_kcal_mol"] == pytest.approx(
        0.20223406721910991, abs=1e-12
    )
    assert metrics["external_validation"]["endpoint_disjoint"]["n"] == 220
    assert metrics["external_validation"]["endpoint_disjoint"]["full_solvai"][
        "mae_kcal_mol"
    ] == pytest.approx(1.152554228782572, abs=1e-12)
    assert metrics["external_validation"]["strict_response_source_disjoint"]["n"] == 97
    assert not table.empty


def test_tier_a_qualification_is_endpoint_disjoint() -> None:
    audit = json.loads(
        (
            ROOT
            / "results/tier_a_external/qualification/tier_a_qualification_summary.json"
        ).read_text()
    )
    assert audit["original_cohort_rows"] == 221
    assert audit["endpoint_overlap_rows"] == 0
    assert audit["endpoint_disjoint_eligible_rows"] == 220
    assert audit["strict_response_source_disjoint_rows"] == 97


def test_manifest_smoke_values_match_test() -> None:
    manifest = json.loads((ROOT / "models/final/manifest.json").read_text())
    card = json.loads((ROOT / "models/final/model_card.json").read_text())
    assert [row["smiles"] for row in manifest["smoke_test"]] == ["CCO", "c1ccccc1"]
    assert card["inference_simulation"] is False
    assert card["input_canonicalization"] == "RDKit canonical isomeric SMILES"
    assert card["calibrated_applicability_score"] is False
    assert card["tier_a_endpoint_disjoint_n"] == 220
    assert card["tier_a_strict_source_disjoint_n"] == 97
    assert card["standardized_teacher_exclusions"] == {
        "combisolv_qm": 2,
        "molsolv_smd": 32,
        "confsolv": 22,
    }
