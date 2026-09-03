# Active SolvAI v3 final gate decision

## Decision: do not scale

Active SolvAI v3 stops before new simulation. The frozen zero-new-simulation
gate found no molecule-held-out evidence that the fifteen SolvAI response
coordinates identify window-level sampling difficulty beyond target-trajectory
diagnostics. A separately frozen power/reference calculation found that the
first adequate predeclared resolving design would require approximately 69,138
reserved RTX 3090 GPU-hours under optimistic noise scaling. This combination is
a quantitative no-go under the blueprint.

## What was frozen before each result

1. `release/V3_MASTER_PROTOCOL_FREEZE.md`, commit `c19101b`, fixed the causal
   contrast, data exposure, target, folds, diagnostics, model ladder, metrics,
   controls, launch criteria and cost accounting before Gate-1 scoring.
2. `release/V3_MASTER_PROTOCOL_FREEZE_AMENDMENT_001.md`, commit `4608c21`,
   clarified the deterministic lambda interactions required for a
   molecule-conditioned linear model before any fit or held-out score.
3. `release/V3_POWER_RELIABILITY_PROTOCOL.md`, commit `9610813`, fixed the
   effect alternatives, inherited heterogeneity, reference model, design grid,
   Monte Carlo procedure, cost rule and adequacy thresholds before power was
   calculated.

## What was actually run

- Complete reproduction of parent SolvAI, Active v1 and Active v2 canonical
  machine-readable results.
- Workspace-wide inventory of 418 `SYSTEM` energy files and all inherited
  scalar response artifacts.
- Deterministic construction of 9,000 post-initial response frames from the
  inherited 12 x 15 x 5-ps dense library.
- Frozen SolvAI inference for the exact fifteen response coordinates on the
  twelve development molecules.
- Nested leave-one-molecule-out Gate-1 analysis over four prefix lengths, three
  fixed target stabilizers, six model variants and five response shuffles.
- 100,000-resample molecule-clustered paired intervals for Gate 1.
- 50,000-trial simulation-based power and reference-reliability calculations
  over the predeclared design grid.
- Deterministic sequential-replay implementation and synthetic/future-access
  tests.

**New v3 production physics:** 0 molecules, 0 windows, 0 trajectories, 0 ps and
0 GPU-hours. No simulation failed because none was scientifically authorized.

## Primary results

| Result | Value |
|---|---:|
| Generic observed diagnostic MAE | 1.427 log-difficulty units |
| SolvAI-conditioned diagnostic MAE | 1.475 |
| SolvAI minus generic | +0.0477 |
| Molecule-clustered 90% CI | -0.0171 to +0.1199 |
| Molecules improved | 6/12 |
| SolvAI relative gain | -3.34% |
| Aligned minus mean shuffled | -0.0102 (90% CI -0.0478 to +0.0310) |
| SolvAI within-molecule rank correlation | 0.225 |

Generic diagnostics reduced error by 7.7% versus lambda/protocol, with a paired
difference of -0.1197 (90% CI -0.1994 to -0.0410). This did not clear the
frozen 10% materiality gate and is not an allocation-policy result.

The SolvAI increment was unfavorable at every prefix: +0.028, +0.048, +0.008
and +0.016 log-difficulty units at 0.5, 1, 2 and 3 ps. It was also unfavorable
under both target-stabilizer sensitivities.

## Independent-reference reliability and minimum resolving design

No compatible inherited molecule-window has independent replicas, so empirical
independent-reference reliability was unavailable. The inherited 2.5-ps
complementary-block dense-integral disagreement had previously been 1.866 kcal
mol-1 on average. Under the frozen optimistic inverse-square-root extrapolation,
the first predeclared design satisfying both policy-effect power and reference
reliability was:

- 256 molecules;
- 15 lambda windows;
- two independent streams;
- 3,000 ps per stream-window;
- 23,040,000 ps total production;
- 62,286.8 projected production GPU-hours;
- 69,138.3 reserved GPU-hours with contingency;
- powers 0.99974 and 0.82970 at the five- and seven-window inherited effect
  scales;
- reference-gate probability 0.82704.

This is the smallest point in the frozen grid and a lower-bound planning result,
not a validated runtime estimate or launch proposal.

## Required handoff answers

1. **What was frozen?** The master protocol, its pre-result feature-map
   clarification and the post-Gate-1 power/reference protocol, all committed
   before their governed results.
2. **What was run?** The inherited reproduction/inventory, Gate-1
   molecule-held-out analysis, destructive controls, deterministic replay tests
   and power/reference planning. No new physics.
3. **Physical work?** Zero new molecules, windows, trajectories, picoseconds or
   GPU-hours.
4. **Did generic adaptation beat equal time?** Not tested as an allocation
   policy because no independent streams were available. Generic diagnostics
   beat lambda/protocol prediction by 7.7%, below the materiality gate.
5. **Did SolvAI-conditioned allocation beat generic?** No policy was authorized;
   its required identifiability precursor was 3.34% worse than generic.
6. **Paired effect and interval?** +0.0477 log-difficulty units, 90% CI -0.0171
   to +0.1199; six of twelve molecules improved.
7. **Independent-reference reliability?** Not empirically available. The first
   planning design with at least 0.80 probability of passing used 3 ns per
   stream-window and 256 molecules.
8. **Coverage and stopping?** Not evaluated because no policy survived to an
   independent replay. No calibration claim is made.
9. **Failures/exclusions?** All twelve dense molecules and 180 windows were
   retained. One unrelated incomplete historical probe remains recorded but was
   not scientifically eligible. There were no v3 simulation failures.
10. **Supported claim?** Early variance diagnostics contain modest development
    signal, but the frozen SolvAI response coordinates do not add held-out
    sampling-difficulty information on the inherited dense library.
11. **Is scale-up justified?** No. The AI gate failed and the minimum resolving
    design is computationally disproportionate.
12. **Next action if the question is reopened?** Treat it as a new program that
    develops simulation-difficulty-specific molecular targets on a much larger
    independent trajectory corpus; do not tune or extend this frozen branch.

## Paths

- Main report: `active_solvai_v3/reports/FINAL_GATE_DECISION.md`
- Evidence audit: `active_solvai_v3/reports/INHERITED_EVIDENCE_AUDIT.md`
- Gate-1 report: `active_solvai_v3/reports/GATE1_IDENTIFIABILITY.md`
- Power report: `active_solvai_v3/reports/POWER_REFERENCE_RELIABILITY_GATE.md`
- Master freeze: `active_solvai_v3/release/V3_MASTER_PROTOCOL_FREEZE.md`
- Main paper: `active_solvai_v3/paper/main.pdf`
- Supplement: `active_solvai_v3/paper/supplementary/supplementary.pdf`
- Canonical metrics: `active_solvai_v3/results/canonical_metrics.json`
- Machine-readable Gate-1 results: `active_solvai_v3/results/gate1/`
- Machine-readable power results: `active_solvai_v3/results/power_reliability/`
- Reproduction: `make -C active_solvai_v3 verify`
