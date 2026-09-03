"""Append-only experiment and compute ledger helpers."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state(root: Path) -> tuple[str, bool]:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = bool(
        subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip()
    )
    return commit, dirty


def append_record(path: Path, record: dict[str, Any]) -> None:
    """Append a single normalized JSON object to a ledger.

    Existing bytes are never rewritten. Callers are responsible for stable run IDs.
    """
    row = dict(record)
    row.setdefault("timestamp_utc", datetime.now(UTC).isoformat())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def assert_unique_run_ids(path: Path) -> None:
    seen: set[str] = set()
    if not path.exists():
        return
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        run_id = json.loads(line)["run_id"]
        if run_id in seen:
            raise AssertionError(f"Duplicate run_id {run_id!r} at ledger line {number}")
        seen.add(run_id)
