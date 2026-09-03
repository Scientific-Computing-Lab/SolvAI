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

**Closed at the quantitative pre-simulation gate.** On the preregistered
12-molecule identifiability test, the matched generic diagnostic model achieved
a log-difficulty MAE of 1.427 and the SolvAI-conditioned model achieved 1.475.
The SolvAI-minus-generic paired difference was +0.0477 (90% CI -0.0171 to
+0.1199), and aligned priors were indistinguishable from five molecule-shuffled
controls. No v3 molecular simulation was launched.

The first numerically adequate point in the separately frozen power/reference
grid requires 256 molecules, 15 windows, two independent streams and 3,000 ps
per stream-window: 23.04 microseconds of production and approximately 69,138
reserved RTX 3090 GPU-hours under optimistic inverse-square-root noise scaling.
This is a quantitative lower-bound design, not a launch recommendation.

## Reproduction

Commands will be added as each prospective gate is frozen. The inherited evidence
reproduction is:

```bash
uv run --project active_solvai python active_solvai_v3/scripts/reproduce_inherited.py
uv run --project active_solvai python active_solvai_v3/scripts/inventory_trajectories.py
PYTHONPATH=active_solvai_v3/src uv run --project . \
  python active_solvai_v3/scripts/build_gate1_dataset.py
PYTHONPATH=active_solvai_v3/src uv run --project active_solvai_v3 \
  python active_solvai_v3/scripts/run_gate1_identifiability.py
PYTHONPATH=active_solvai_v3/src uv run --project active_solvai_v3 \
  python active_solvai_v3/scripts/run_power_reliability_gate.py
PYTHONPATH=active_solvai_v3/src uv run --project active_solvai_v3 \
  python active_solvai_v3/scripts/finalize_v3.py
make -C active_solvai_v3 test paper
```

The main decision report is
[`reports/FINAL_GATE_DECISION.md`](reports/FINAL_GATE_DECISION.md).
