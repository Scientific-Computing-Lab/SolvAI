#!/usr/bin/env python3
"""Verify model hashes and reproduce the frozen SMILES-only smoke predictions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from solv_ai import predict_smiles

ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = ROOT / "models/final"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((MODEL_DIR / "manifest.json").read_text())
    for filename, expected in manifest["model_files"].items():
        path = MODEL_DIR / filename
        if not path.is_file() or path.stat().st_size != expected["bytes"]:
            raise AssertionError(f"Missing or truncated artifact: {filename}")
        observed = sha256(path)
        if observed != expected["sha256"]:
            raise AssertionError(f"SHA-256 mismatch: {filename}")
    smoke = manifest["smoke_test"]
    predictions, spread = predict_smiles([row["smiles"] for row in smoke])
    expected_prediction = np.asarray([row["prediction_kcal_mol"] for row in smoke])
    expected_spread = np.asarray([row["ensemble_spread_kcal_mol"] for row in smoke])
    if not np.allclose(predictions, expected_prediction, atol=2e-6, rtol=0):
        raise AssertionError("Frozen predictions were not reproduced")
    if not np.allclose(spread, expected_spread, atol=2e-6, rtol=0):
        raise AssertionError("Frozen ensemble spreads were not reproduced")
    result = {
        "status": "PASS",
        "files_verified": len(manifest["model_files"]),
        "smiles_only_predictions_verified": len(smoke),
        "simulation_at_inference": False,
        "tolerance": 2e-6,
        "input": "SMILES only",
        "forbidden_inputs": [
            "experimental hydration free energy",
            "ARROW or PIMD free energy",
            "trajectory or probe output",
            "family, scaffold or fold metadata",
        ],
    }
    audit_dir = ROOT / "audits"
    audit_dir.mkdir(parents=True, exist_ok=True)
    (audit_dir / "artifact_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    rows = [
        {"artifact": filename, **values, "status": "PASS"}
        for filename, values in manifest["model_files"].items()
    ]
    pd.DataFrame(rows).to_csv(audit_dir / "artifact_audit.csv", index=False)
    lines = [
        "# SolvAI inference-artifact audit",
        "",
        "**Status: PASS.** All packaged file hashes and two end-to-end SMILES-only",
        "predictions reproduce within 2 × 10⁻⁶ kcal/mol.",
        "",
        f"- Files verified: {len(rows)}",
        "- Input: SMILES only",
        "- Simulation at inference: no",
        "- Benchmark or prediction table read by runtime: no",
        "- Experimental, ARROW, PIMD, trajectory, probe, family, scaffold or fold input: no",
    ]
    (audit_dir / "artifact_audit.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
