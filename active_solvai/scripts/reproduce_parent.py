#!/usr/bin/env python3
"""Reproduce the minimum frozen SolvAI evidence without rewriting parent outputs."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai/results/phase0/parent_reproduction.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((ROOT / "models/final/manifest.json").read_text())
    hashes = {}
    for name, expected in manifest["model_files"].items():
        path = ROOT / "models/final" / name
        observed = sha256(path)
        if observed != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            raise AssertionError(f"Parent artifact mismatch: {name}")
        hashes[name] = observed

    smoke = manifest["smoke_test"]
    # Use the parent's locked runtime because Active SolvAI intentionally has a
    # smaller analysis environment. The subprocess reads the artifact only.
    parent_python = ROOT / ".venv/bin/python"
    snippet = (
        "import json; from solv_ai import predict_smiles; "
        f"p,s=predict_smiles({[row['smiles'] for row in smoke]!r}); "
        "print(json.dumps({'prediction':p.tolist(),'spread':s.tolist()}))"
    )
    environment = os.environ.copy()
    environment["PATH"] = f"{ROOT / '.venv/bin'}:{environment.get('PATH', '')}"
    runtime = json.loads(
        subprocess.check_output(
            [str(parent_python), "-c", snippet], cwd=ROOT, text=True, env=environment
        ).strip()
    )
    prediction = np.asarray(runtime["prediction"], dtype=float)
    spread = np.asarray(runtime["spread"], dtype=float)
    expected_prediction = np.array([row["prediction_kcal_mol"] for row in smoke])
    expected_spread = np.array([row["ensemble_spread_kcal_mol"] for row in smoke])
    if not np.allclose(prediction, expected_prediction, atol=2e-6, rtol=0):
        raise AssertionError("Parent smoke predictions did not reproduce")
    if not np.allclose(spread, expected_spread, atol=2e-6, rtol=0):
        raise AssertionError("Parent smoke spreads did not reproduce")

    frame = pd.read_parquet(
        ROOT / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    primary = frame.loc[frame["partition"].eq("standardized_exclusion_primary")]
    maes = primary.groupby("method", observed=True)["absolute_error"].mean().to_dict()
    expected_metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    expected_structure = float(expected_metrics["methods"]["matched_structure_only"]["mae_kcal_mol"])
    expected_full = float(expected_metrics["methods"]["full_solvai"]["mae_kcal_mol"])
    if not np.isclose(maes["A_structure_only"], expected_structure, atol=1e-12):
        raise AssertionError("Parent structure-only metric did not reproduce")
    if not np.isclose(maes["F_full_solvai"], expected_full, atol=1e-12):
        raise AssertionError("Parent SolvAI metric did not reproduce")

    payload = {
        "status": "PASS",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "parent_commit_expected": "531f6cfd21e319c951b461c9ef24fa754790f91d",
        "model_files_verified": len(hashes),
        "model_hashes": hashes,
        "smoke_predictions": [
            {"smiles": row["smiles"], "prediction_kcal_mol": float(p), "spread": float(s)}
            for row, p, s in zip(smoke, prediction, spread, strict=True)
        ],
        "fixed_primary": {
            "n": int(primary["molecule_id"].nunique()),
            "structure_only_mae_kcal_mol": float(maes["A_structure_only"]),
            "solvai_mae_kcal_mol": float(maes["F_full_solvai"]),
        },
        "source_hashes": {
            "paper_metrics_json": sha256(ROOT / "results/paper_metrics.json"),
            "parent_predictions": sha256(
                ROOT / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
            ),
            "main_pdf": sha256(ROOT / "paper/main.pdf"),
            "supplement_pdf": sha256(ROOT / "paper/supplementary/supplementary.pdf"),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["fixed_primary"], indent=2))


if __name__ == "__main__":
    main()
