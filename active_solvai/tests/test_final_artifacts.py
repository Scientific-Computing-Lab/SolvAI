from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_final_decision_is_consistent_with_frozen_primary_comparisons() -> None:
    result = json.loads((ROOT / "results/phase2/active_solvai_final_metrics.json").read_text())
    assert result["program_decision"] == "NO_GO"
    assert result["phase1_endpoint_decision"] == "negative"
    assert result["phase2_decision"] == "negative"
    assert result["direction_status"]["C_adaptive_multifidelity"].startswith("not launched")
    for row in result["phase2_primary_comparisons"]:
        assert row["budget"] in (5, 7)
        assert row["active_mae"] > row["comparator_mae"]
        assert not row["condition_one_passed"]


def test_prospective_compute_accounting_is_complete() -> None:
    result = json.loads((ROOT / "results/phase2/active_solvai_final_metrics.json").read_text())
    prospective = result["simulation_compute"]["by_role"]["prospective"]
    assert prospective["attempted_windows"] == 96
    assert prospective["passed_windows"] == 96
    assert prospective["failed_windows"] == 0
    assert prospective["production_ps"] == 480
