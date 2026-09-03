# Active SolvAI

Active SolvAI is an isolated follow-up to the completed SolvAI release. It asks whether a small number of **actual target-molecule solvation-response observations** can update the frozen structure-only SolvAI prior and reduce the cost of reconstructing an alchemical response or predicting experimental hydration free energy.

The parent SolvAI release is immutable. Active SolvAI lives on the `active-solvai` branch and writes only below this directory.

## Scientific targets

1. **Same-Hamiltonian reconstruction:** infer hidden alchemical response windows and the dense integral under the simulation protocol that generated the observations.
2. **Experimental endpoint:** test whether actual response information improves held-out experimental hydration free energies beyond frozen SolvAI.

These targets are evaluated and reported separately.

## Gated workflow

```text
Phase 0: reproduce parent + inventory existing response data
  -> Phase 1: nested actual-observation gate (endpoint route negative;
              sparse reconstruction diagnostic conditionally positive)
  -> Phase 2: bounded prospective dense sentinel for reconstruction
  -> Phase 3: prospective sentinel pilot, only after replay is frozen and passes
  -> Tier-B: only after an explicit scale request and PI approval
```

The authoritative specification is `Active_SolvAI_Master_Research_Blueprint.docx` at the workspace root. Its ingestion record is in [`reports/ACTIVE_SOLVAI_BLUEPRINT_INGESTION.md`](reports/ACTIVE_SOLVAI_BLUEPRINT_INGESTION.md).

## Reproduction

```bash
uv sync --project active_solvai --locked
uv run --project active_solvai pytest active_solvai/tests -q
uv run --project active_solvai python active_solvai/scripts/reproduce_parent.py
uv run --project active_solvai python active_solvai/scripts/inventory_responses.py
uv run --project active_solvai python active_solvai/scripts/run_phase1_gate.py
uv run --project active_solvai python active_solvai/scripts/summarize_phase1.py
```

The Phase 1 experimental-endpoint gate is a registered null result: actual
three-point PIMD2 residuals increased five-repeat MAE by 0.00348 kcal mol⁻¹
relative to frozen SolvAI. A separate held-point reconstruction diagnostic
supports one bounded dense same-Hamiltonian sentinel test; it does not reopen
the endpoint claim.

Decisive analyses may run only after the corresponding freeze has been committed. Every attempted run is appended to `runs/ledger.jsonl`; failed work is retained and counted.
