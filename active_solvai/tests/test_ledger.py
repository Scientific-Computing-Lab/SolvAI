from __future__ import annotations

import json

import pytest

from active_solvai.ledger import append_record, assert_unique_run_ids


def test_append_only_ledger_and_duplicate_detection(tmp_path):
    path = tmp_path / "ledger.jsonl"
    append_record(path, {"run_id": "r1", "status": "started"})
    append_record(path, {"run_id": "r2", "status": "complete"})
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    assert [row["run_id"] for row in rows] == ["r1", "r2"]
    assert_unique_run_ids(path)
    append_record(path, {"run_id": "r1", "status": "retry"})
    with pytest.raises(AssertionError, match="Duplicate run_id"):
        assert_unique_run_ids(path)
