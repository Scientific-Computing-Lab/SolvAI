# Active SolvAI v3: Adaptive Simulation Effort

This isolated workstream tests whether frozen, molecule-conditioned SolvAI response
coordinates improve allocation of **additional trajectory time** across the fixed
15-window ARROW/PIMD2 alchemical grid. The primary comparison is a SolvAI-conditioned
allocator versus a matched generic allocator that sees only the target simulation's
observed data. Equal-time and other predeclared controls are required.

The completed SolvAI paper and Active SolvAI v1/v2 no-go branches are immutable.
Experimental hydration labels are not policy inputs. Development decisions are scored
only against complementary or independent trajectory streams, and a prospective
sentinel may begin only after its policy and cohort freezes are committed.

## Status

Phase 0 is in progress. No v3 molecular simulation has been launched.

## Reproduction

Commands will be added as each prospective gate is frozen. The inherited evidence
reproduction is:

```bash
uv run --project active_solvai python active_solvai_v3/scripts/reproduce_inherited.py
uv run --project active_solvai python active_solvai_v3/scripts/inventory_trajectories.py
```

