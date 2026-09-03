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

## Ledger rules

- Append a record before and after every attempted analysis or simulation run.
- Never aggregate by discovering filenames; aggregation reads `runs/ledger.jsonl`.
- Failed and aborted work remains in the ledger and compute totals.
- A decisive scoring command must refuse to run unless its freeze path exists and is committed.

Phase-specific commands will be added before their corresponding freeze is committed.

