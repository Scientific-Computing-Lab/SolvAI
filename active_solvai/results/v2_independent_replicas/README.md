# Independent-replica resolution outputs

This directory records the prospectively gated attempt to resolve whether stable molecule-specific lambda placement exists beyond short-trajectory noise.

The gate failed before simulation. Consequently there are no new trajectories or molecule-level independent-replica predictions. `resolution_status.json` records this explicitly; the complete Monte Carlo outputs are under `power/`.

Reproduce the gate with:

```bash
uv run --project active_solvai python active_solvai/scripts/run_independent_replica_power.py
```
