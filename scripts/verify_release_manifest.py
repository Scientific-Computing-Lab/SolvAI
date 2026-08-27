#!/usr/bin/env python3
"""Verify every canonical artifact frozen in the release manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    manifest = json.loads((ROOT / "release_manifest.json").read_text())
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file() or path.stat().st_size != expected["bytes"]:
            raise AssertionError(f"Missing or truncated release artifact: {relative}")
        if sha256(path) != expected["sha256"]:
            raise AssertionError(f"Release artifact hash mismatch: {relative}")
    print(f"Release manifest PASS ({len(manifest['files'])} canonical artifacts)")


if __name__ == "__main__":
    main()
