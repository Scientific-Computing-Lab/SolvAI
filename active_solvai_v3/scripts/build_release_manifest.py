#!/usr/bin/env python3
"""Hash the final v3 release artifacts after all deterministic builds."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ACTIVE = ROOT / "active_solvai_v3"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    patterns = (
        "release/*.md",
        "reports/*.md",
        "results/**/*.json",
        "results/**/*.csv",
        "results/**/*.parquet",
        "figures/*.pdf",
        "figures/*.svg",
        "figures/*.png",
        "paper/main.tex",
        "paper/main.pdf",
        "paper/generated_metrics.tex",
        "paper/supplementary/supplementary.tex",
        "paper/supplementary/supplementary.pdf",
        "runs/ledger.jsonl",
    )
    paths = sorted({path for pattern in patterns for path in ACTIVE.glob(pattern) if path.is_file()})
    payload = {
        "schema_version": 1,
        "branch": subprocess.check_output(
            ["git", "branch", "--show-current"], cwd=ROOT, text=True
        ).strip(),
        "source_parent_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "artifacts": [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in paths
        ],
    }
    output = ACTIVE / "release/V3_ARTIFACT_MANIFEST.json"
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {output} with {len(paths)} artifacts")


if __name__ == "__main__":
    main()
