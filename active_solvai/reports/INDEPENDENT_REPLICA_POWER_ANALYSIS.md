# Power and precision gate for independent PIMD2 replicas

**Prospective protocol:** commit `52827f851ed0d1a540b0492c63c3df269cf774e1`  
**Decision:** FAIL — the registered eight-molecule simulation must not launch.

This calculation precedes every new simulation. It asks whether two 50-ps replicas on the eight frozen sentinels can resolve a true cross-replica oracle advantage of 0.30 kcal mol-1 at both five and seven windows while also producing a reliable dense reference.

## Registered eight-molecule design

| Quantity | Result |
|---|---:|
| Power, 5 windows | 0.000 |
| Power, 7 windows | 0.000 |
| Conservative power, 5 windows | 0.000 |
| Conservative power, 7 windows | 0.001 |
| Probability dense-reference gate passes | 0.446 |
| Required probability for each primary gate | 0.800 |
| Projected production GPU-hours | 32.44 |
| Operational reservation GPU-hours | 36.01 |

The design is launchable only when both budget-specific powers and dense reliability are at least 0.80. The decision is therefore mechanical rather than result-dependent.

## Smallest adequately powered alternative

No design on the frozen grid passed all power and reliability criteria.

This alternative is **not authorized or launched**. If it requires more than eight molecules, those molecules must be chosen under a separate chemistry-first prospective freeze; this power calculation does not select them.

No amount of additional trajectory length or molecule replication on the registered empirical effect distribution makes the *full* decision rule adequately powered: after centering to a mean benefit of 0.30 kcal mol-1, only 4/8 molecules at five windows and 5/8 at seven windows retain favorable effects, below the required 75% consistency. Increased sample size therefore estimates that inconsistency more precisely rather than curing it.

For planning context only, dropping the consistency requirement and asking solely whether the **mean** effect is below zero gives a first passing boundary at **105 molecules, two replicas, 15 windows and 75 ps per replica**. Its simulated powers are 0.801 and 0.817, with dense-reference reliability probability 0.940. It requires 638.7 projected production GPU-h and 708.9 reserved GPU-h. This would resolve an average effect only; it cannot establish stable, broadly shared molecule-specific placement and is not recommended or authorized.

## Interpretation

The analysis preserves the observed between-molecule effect pattern, shifts only its mean to the scientifically meaningful alternative of -0.30 kcal mol-1, and projects finite-trajectory residual variance from 2.5 ps to the candidate duration by inverse-time scaling. The conservative sensitivity doubles residual variance. Dense-integral reliability is independently projected from the observed complementary-half discrepancies.

The calculation cannot prove that inverse-time scaling will hold or that a true oracle benefit exists. It only determines whether the proposed experiment has a reasonable chance to resolve an effect worth pursuing under the most favorable registered noise law. Failure therefore blocks simulation; passing merely authorizes the frozen experiment.

## Reproduction

```bash
uv run --project active_solvai python active_solvai/scripts/run_independent_replica_power.py
```

Machine-readable results are in `active_solvai/results/v2_independent_replicas/power/`.
