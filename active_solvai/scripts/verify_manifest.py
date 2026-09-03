#!/usr/bin/env python3
"""Verify every file recorded in the Active SolvAI SHA-256 manifest."""

from __future__ import annotations

import json
from pathlib import Path

from active_solvai.ledger import sha256

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.json"


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    failures: list[str] = []
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        observed_hash = sha256(path)
        observed_size = path.stat().st_size
        if observed_hash != expected["sha256"]:
            failures.append(f"hash mismatch: {relative}")
        if observed_size != expected["bytes"]:
            failures.append(f"size mismatch: {relative}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"verified {len(manifest['files'])} files from {MANIFEST}")


if __name__ == "__main__":
    main()
