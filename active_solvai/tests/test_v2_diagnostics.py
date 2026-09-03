from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/v2_diagnostics/oracle_independent_noise"


def test_cross_block_design_is_complete_and_disjoint() -> None:
    predictions = pd.read_parquet(RESULTS / "cross_block_predictions.parquet")
    assert len(predictions) == 8 * 6 * 7 * 2
    assert set(predictions.total_windows) == {5, 7}
    assert predictions.molecule_id.nunique() == 8
    assert predictions.method.nunique() == 7
    for selection, evaluation in (
        predictions[["selection_blocks", "evaluation_blocks"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ):
        left = {int(value) for value in selection.split(",")}
        right = {int(value) for value in evaluation.split(",")}
        assert left.isdisjoint(right)
        assert left | right == {0, 1, 2, 3}


def test_registered_diagnostic_conclusion() -> None:
    metrics = json.loads((RESULTS / "canonical_metrics.json").read_text())
    assert metrics["conclusion"] == "C. INCONCLUSIVE DUE TO RESPONSE NOISE"
    assert metrics["reliability_gate_exceeded"] is True
    assert metrics["b_rule_passed"] is False
    assert metrics["a_rule_passed"] is False
