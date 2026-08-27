#!/usr/bin/env python3
"""Create the immutable checksum and smoke-prediction manifest for the release model."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

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
    model_files = sorted(
        path for path in MODEL_DIR.iterdir() if path.is_file() and path.name != "manifest.json"
    )
    smoke_smiles = ["CCO", "c1ccccc1"]
    prediction, spread = predict_smiles(smoke_smiles)
    manifest = {
        "schema_version": 1,
        "model_files": {
            path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in model_files
        },
        "smoke_test": [
            {
                "smiles": smiles,
                "prediction_kcal_mol": float(value),
                "ensemble_spread_kcal_mol": float(uncertainty),
            }
            for smiles, value, uncertainty in zip(smoke_smiles, prediction, spread, strict=True)
        ],
    }
    (MODEL_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
