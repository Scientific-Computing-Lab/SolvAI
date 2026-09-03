import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def _load_module():
    path = ROOT / "active_solvai/scripts/run_independent_replica_power.py"
    spec = importlib.util.spec_from_file_location("independent_replica_power", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_power_inputs_and_centering() -> None:
    module = _load_module()
    assert module.assert_inputs() == {key: expected for key, (_, expected) in module.INPUTS.items()}
    patterns, residual_variances, components = module.effect_components()
    assert set(patterns) == {5, 7}
    assert all(len(values) == 8 for values in patterns.values())
    assert all(np.isclose(values.mean(), 0.0) for values in patterns.values())
    assert all(value > 0 for value in residual_variances.values())
    assert len(components) == 16


def test_frozen_gate_blocks_simulation() -> None:
    result = json.loads(
        (
            ROOT
            / "active_solvai/results/v2_independent_replicas/power/power_analysis.json"
        ).read_text()
    )
    assert result["protocol_commit"] == "52827f851ed0d1a540b0492c63c3df269cf774e1"
    assert result["launch_authorized_by_gate"] is False
    assert result["simulation_started"] is False
    assert result["proposed_design"]["power_budget_5"] < 0.80
    assert result["proposed_design"]["power_budget_7"] < 0.80
    assert result["proposed_design"]["dense_reliability_probability"] < 0.80
