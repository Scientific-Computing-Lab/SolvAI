#!/usr/bin/env python3
"""Build a deterministic SHA-256 manifest for tracked Active SolvAI artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "MANIFEST.json"
EXCLUDED_PARTS = {
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    "FFDbLog",
    "cache",
    "output",
}
EXCLUDED_SUFFIXES = {
    ".aux",
    ".bbl",
    ".bcf",
    ".blg",
    ".fdb_latexmk",
    ".fls",
    ".log",
    ".out",
    ".run.xml",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    files = {}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == OUT:
            continue
        relative = path.relative_to(ROOT)
        if EXCLUDED_PARTS.intersection(relative.parts):
            continue
        if any(path.name.endswith(suffix) for suffix in EXCLUDED_SUFFIXES):
            continue
        files[str(relative)] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    OUT.write_text(json.dumps({"schema_version": 1, "files": files}, indent=2) + "\n")
    print(f"wrote {OUT} ({len(files)} files)")


if __name__ == "__main__":
    main()
