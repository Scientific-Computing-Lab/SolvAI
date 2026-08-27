#!/usr/bin/env python3
"""Hash the canonical scientific and publication artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = [
    "data/benchmark/arrow_solvation_master.parquet",
    "data/manifests/endpoint_label_manifest.json",
    "data/manifests/training_source_manifest.json",
    "results/paper_metrics.json",
    "results/paper_metrics.csv",
    "results/predictions/headline_oof.parquet",
    "results/predictions/hard_holdout_oof.parquet",
    "results/robustness/repeated_oof.parquet",
    "audits/leakage_audit.json",
    "audits/artifact_audit.json",
    "audits/claim_red_team.json",
    "audits/security_audit.json",
    "models/final/manifest.json",
    "results/runtime/runtime_benchmark.json",
    "paper/main.tex",
    "paper/main.pdf",
    "paper/supplementary.tex",
    "paper/supplementary.pdf",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = {}
    for relative in TARGETS:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        files[relative] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {"schema_version": 1, "files": files}
    output = ROOT / "release_manifest.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Hashed {len(files)} canonical artifacts")


if __name__ == "__main__":
    main()
