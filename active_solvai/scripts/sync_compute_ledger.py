"""Materialize the append-only JSONL run ledger as a reviewable CSV."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "runs/ledger.jsonl"
OUTPUT = ROOT / "runs/COMPUTE_LEDGER.csv"

COLUMNS = [
    "run_id",
    "parent_run_id",
    "timestamp_utc",
    "git_commit",
    "dirty",
    "config_sha256",
    "molecule_id",
    "exposure_tier",
    "outer_split",
    "action",
    "stage",
    "status",
    "fidelity",
    "lambda",
    "beads",
    "requested_steps",
    "completed_steps",
    "seed",
    "device",
    "wall_seconds",
    "wall_hours",
    "gpu_hours",
    "cpu_hours",
    "simulated_time_ps",
    "force_evaluations",
    "bead_windows",
    "qc_status",
    "quality_control",
    "failure_reason",
    "input_sha256",
    "output_sha256",
    "freeze_commit",
    "command",
]


def main() -> None:
    records = [json.loads(line) for line in SOURCE.read_text().splitlines() if line.strip()]
    run_ids = [record["run_id"] for record in records]
    if len(run_ids) != len(set(run_ids)):
        raise AssertionError("Duplicate run IDs in append-only ledger")
    frame = pd.DataFrame(records)
    for column in COLUMNS:
        if column not in frame:
            frame[column] = None
    frame[COLUMNS].to_csv(OUTPUT, index=False)
    print(f"wrote {OUTPUT} with {len(frame)} runs")


if __name__ == "__main__":
    main()
