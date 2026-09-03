#!/usr/bin/env python3
"""Verify frozen v3 results, artifact completeness and claim consistency."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ACTIVE = ROOT / "active_solvai_v3"


def main() -> None:
    canonical = json.loads((ACTIVE / "results/canonical_metrics.json").read_text())
    gate = json.loads((ACTIVE / "results/gate1/gate1_canonical_metrics.json").read_text())
    power = json.loads(
        (ACTIVE / "results/power_reliability/power_reliability_canonical.json").read_text()
    )
    predictions = pd.read_parquet(ACTIVE / "results/gate1/gate1_oof_predictions.parquet")

    assert canonical["decision"] == "quantitative_gate_failed_no_simulation"
    assert canonical["new_simulation"] == {
        "molecules": 0,
        "trajectories": 0,
        "simulated_ps": 0.0,
        "gpu_hours": 0.0,
    }
    assert not gate["generic_gate_passed"]
    assert not gate["solvai_identifiability_gate_passed"]
    assert not power["launch_authorized"]
    assert len(predictions) == 23_760
    key = ["model", "stabilizer", "prefix_ps", "molecule_id", "lambda"]
    assert not predictions.duplicated(key).any()
    assert predictions.molecule_id.nunique() == 12
    assert np.isfinite(predictions[["target", "prediction"]]).all().all()

    required = [
        ACTIVE / "release/V3_MASTER_PROTOCOL_FREEZE.md",
        ACTIVE / "release/V3_MASTER_PROTOCOL_FREEZE_AMENDMENT_001.md",
        ACTIVE / "release/V3_POWER_RELIABILITY_PROTOCOL.md",
        ACTIVE / "reports/INHERITED_EVIDENCE_AUDIT.md",
        ACTIVE / "reports/GATE1_IDENTIFIABILITY.md",
        ACTIVE / "reports/POWER_REFERENCE_RELIABILITY_GATE.md",
        ACTIVE / "reports/FINAL_GATE_DECISION.md",
        ACTIVE / "paper/main.pdf",
        ACTIVE / "paper/supplementary/supplementary.pdf",
        ACTIVE / "figures/gate1_identifiability.pdf",
        ACTIVE / "figures/v3_power_reliability.pdf",
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    assert not missing, f"missing release artifacts: {missing}"

    forbidden_claims = [
        r"SolvAI-conditioned allocation beats",
        r"independent-replica validation was achieved",
        r"scale-up is justified",
    ]
    text_paths = [
        ACTIVE / "paper/main.tex",
        ACTIVE / "paper/supplementary/supplementary.tex",
        ACTIVE / "README.md",
    ]
    text = "\n".join(path.read_text() for path in text_paths)
    for pattern in forbidden_claims:
        assert not re.search(pattern, text, flags=re.IGNORECASE), pattern
    print("Active SolvAI v3 release verification passed")


if __name__ == "__main__":
    main()
