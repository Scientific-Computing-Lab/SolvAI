# Active SolvAI runbook

All commands are run from the parent repository root.

## Environment

```bash
uv sync --project active_solvai --locked
uv run --project active_solvai pytest active_solvai/tests -q
```

## Phase 0

```bash
uv run --project active_solvai python active_solvai/scripts/capture_environment.py
uv run --project active_solvai python active_solvai/scripts/reproduce_parent.py
uv run --project active_solvai python active_solvai/scripts/inventory_responses.py
```

## Phase 1

The scoring protocol was frozen at commit `a0dd986`.

```bash
uv run --project active_solvai python active_solvai/scripts/run_phase1_gate.py
uv run --project active_solvai python active_solvai/scripts/summarize_phase1.py
```

Canonical results are in `results/phase1/phase1_canonical_metrics.json`; the
human-readable interpretation is `reports/PHASE1_ACTUAL_OBSERVATION_GATE.md`.

## Ledger rules

- Append a record before and after every attempted analysis or simulation run.
- Never aggregate by discovering filenames; aggregation reads `runs/ledger.jsonl`.
- Failed and aborted work remains in the ledger and compute totals.
- A decisive scoring command must refuse to run unless its freeze path exists and is committed.

Phase-specific commands will be added before their corresponding freeze is committed.

## Phase 2 dense sentinel

The scientific protocol was frozen at commit `512aad3`; the calibration lock
was committed at `3e69d10` before prospective simulation.

```bash
uv run --project active_solvai python active_solvai/scripts/collect_dense_sentinel.py --role prospective
uv run --project active_solvai python active_solvai/scripts/evaluate_dense_sentinel.py
uv run --project active_solvai python active_solvai/scripts/summarize_dense_sentinel.py
uv run --project active_solvai python active_solvai/scripts/finalize_negative_program.py
uv run --project active_solvai python active_solvai/scripts/sync_compute_ledger.py
uv run --project active_solvai python active_solvai/scripts/build_manifest.py
uv run --project active_solvai python active_solvai/scripts/verify_manifest.py
```

The result is a registered no-go. Direction C and Tier-B must not be launched
from this branch without a new prospective protocol and explicit authorization.
