# V3 power and independent-reference gate

## Decision

**Do not launch the development or sentinel simulation campaign.** The frozen
SolvAI-conditioned difficulty model failed Gate 1, and the conservative
simulation-based calculation shows that a policy-effect experiment resolving
the registered 20% increment would require a scale far beyond the local
exploratory envelope.

The first numerically adequate design in the prospectively frozen grid is:

| Quantity | Design |
|---|---:|
| Molecules | 256 |
| Lambda windows | 15 |
| Independent streams | 2 |
| Production per stream-window | 3,000 ps |
| Total production | 23,040,000 ps (23.04 microseconds) |
| Projected production work | 62,286.8 RTX 3090 GPU-hours |
| Operational reservation (11% contingency) | **69,138.3 GPU-hours** |
| Power, 5-window 20% effect | 0.99974 |
| Power, 7-window 20% effect | 0.82970 |
| Reference-reliability probability | 0.82704 |

This is the smallest **predeclared grid point**, not a recommendation to spend
that compute. It is based on optimistic inverse-square-root noise scaling and
inherits effect heterogeneity from the failed v1 lambda-placement policy because
no independent v3 effort-allocation effect exists yet. The estimate is therefore
a lower-bound planning diagnostic.

## Frozen calculation

The protocol was committed in
`release/V3_POWER_RELIABILITY_PROTOCOL.md` at `9610813` before this calculation.
It used 50,000 Monte Carlo cohorts per condition with seed 20260903.

The minimum-useful effects were fixed at a 20% reduction relative to the
inherited generic-BQ error scale: -0.341905 kcal mol-1 at five windows and
-0.222900 kcal mol-1 at seven windows. Molecule heterogeneity was obtained by
centering the eight prospective v1 Active-minus-generic paired effects. The
null and a 30% larger alternative were also simulated. Detection required the
upper endpoint of a two-sided 95% paired interval to be below zero.

Reference planning resampled the 24 unique complementary 2.5-ps dense-integral
disagreements from v2 and scaled them as `sqrt(2.5/T)`. The master criterion
required the cohort 90th percentile to be no more than half the smaller
minimum-useful effect, 0.111450 kcal mol-1.

## Why a small campaign cannot answer the question

| Candidate | 5-window power | 7-window power | Reference pass probability | Reserved GPU-hours |
|---|---:|---:|---:|---:|
| 12 molecules, 50 ps | 0.218 | 0.114 | 0.000 | 54.0 |
| 20 molecules, 50 ps | 0.326 | 0.142 | 0.000 | 90.0 |
| 20 molecules, 3,000 ps | 0.326 | 0.142 | 0.734 | 5,401.4 |
| 64 molecules, 3,000 ps | 0.769 | 0.317 | 0.716 | 17,284.6 |
| 128 molecules, 3,000 ps | 0.970 | 0.545 | 0.821 | 34,569.2 |
| 256 molecules, 3,000 ps | 1.000 | 0.830 | 0.827 | 69,138.3 |

More duration repairs only the reference. It does not repair the molecule-level
power deficit at the more difficult seven-window effect. More molecules repair
power but do not by themselves stabilize a short dense reference. Both axes
must increase.

## Scientific conclusion

The present evidence supports a narrow generic observation: variance/SEM-like
prefix diagnostics contain modest information about later sampling difficulty.
It does **not** support the specific Active SolvAI claim that the frozen response
coordinates improve effort allocation beyond a matched generic allocator.

Accordingly:

- no `V3_POLICY_FREEZE.md` is created because no SolvAI-conditioned candidate
  survived the registered identifiability gate;
- no `V3_SENTINEL_FREEZE.md` is created because a prospective sentinel is not
  scientifically authorized;
- no new trajectory, PIMD4/PIMD8 calculation or experimental-endpoint analysis
  is launched;
- the exact resolving scale is documented rather than replaced with an
  underpowered token experiment.

## Reproduction

```bash
PYTHONPATH=active_solvai_v3/src uv run --project active_solvai_v3 \
  python active_solvai_v3/scripts/run_power_reliability_gate.py
```

Machine-readable outputs are in `results/power_reliability/`; the figure source
is the script above and the rendered figure is `figures/v3_power_reliability.*`.
