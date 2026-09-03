# Prospective dense PIMD2 sentinel freeze

**Frozen:** 2026-09-03 UTC  
**Starting scientific commit:** `c84ca60617074d44bafb622a13913ec83042a18f`  
**Campaign:** `AS-P2-DENSE-PIMD2-001`  
**Configuration:** `active_solvai/configs/dense_sentinel_v1.json`  
**Configuration SHA-256:** `fc6afd2a4f8a255a641228df624769e6ed437dc4cc69dc244480a33c2936594d`  
**Prepared manifest:** `active_solvai/simulations/dense_pimd2/manifest.csv`  
**Prepared-manifest SHA-256:** `9c46ad5b78af027f3c82d84ff300c1ea839fbd13650eba704bbfef95f5871358`

This file was committed before any of the 144 previously unavailable dense-window
responses were generated or read. The existing 0.1, 0.5 and 0.9 observations and
all Phase 1 outcomes were already known. This is a prospective test only with
respect to the twelve additional lambda windows per molecule.

## Question and permissible claim

The sole primary question is whether a molecule-conditioned response prior can
reconstruct a **dense, short-window PIMD2 alchemical response and its numerical
integral** from fewer same-Hamiltonian windows than fixed, random and generic
Bayesian alternatives. The experimental-endpoint route was killed in Phase 1 and
will not be reopened with these molecules.

Success can support a claim about short-window, same-Hamiltonian response
reconstruction. It cannot establish convergence to full PIMD8, agreement with
experiment, chemical blindness, unbiased thermodynamic integration, or a
production speedup.

## Molecules fixed before dense responses

Selection used chemistry, prior parameterization success and breadth only. It did
not use per-molecule SolvAI errors or ungenerated dense responses.

### Calibration set (new responses may tune only frozen calibration choices)

| Molecule | Family | Role |
|---|---|---|
| Propane | alkane | small non-polar |
| Ethanol | alcohol | hydrogen-bond donor/acceptor |
| Acetone | ketone | polar acceptor |
| pyridine | heteroaromatic | aromatic acceptor |

### Prospective sentinel set (eight-molecule primary evaluation)

| Molecule | Family | Diagnostic role |
|---|---|---|
| Cyclohexane | alkane | cyclic non-polar |
| Octane | alkane | flexible hydrophobe |
| Octanol | alcohol | large amphiphile |
| AceticAcid | acid | strong hydrogen bonding |
| Acetamide | amide | difficult polar neutral |
| DiMethylEther | ether | compact acceptor |
| EthylAcetate | ester | flexible polar neutral |
| Anthracene | aromatic | fused aromatic |

All twelve are ARROW-development molecules with already viewed experimental
labels, three-point PIMD2 responses and parent-model predictions. The primary
evaluation is therefore a prospective simulation-response sentinel, not an
external or blind molecular cohort.

## Simulation protocol

- Native Arbalest executable:
  `/home/galoren/Freecurve_AI_Solvation/repositories/arbalest/build_solvation_gcc10/ARBALEST/ARBALEST`.
- ARROW force field `933`; explicit water; 24 Å periodic cube; 298 K and 1 bar.
- Native PIMD2 with PILE thermostat and the inherited Berendsen production
  barostat; 0.002 ps time step.
- Per window: 250 minimization steps, 500 equilibration steps and 2,500
  production steps (5 ps).
- Fifteen fixed lambda values:
  `0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0`.
- Existing matched observations at 0.1, 0.5 and 0.9 are reused. Twelve new
  windows are generated for each molecule. The full common pool is 180 windows,
  of which 144 are new.
- Production energy output is retained every 50 steps. Coordinate trajectories
  are disabled only for the new windows to reduce storage; this does not change
  forces, integration or energy output.
- A new window passes quality control only if Arbalest exits zero, prints its
  success marker, produces exactly one finite SYSTEM energy series with at least
  45 frames, and its mean temperature and density lie in [250, 350] K and
  [0.70, 1.30] g cm^-3. Failures remain in the ledger and are rerun at most once
  without changing the configuration. A second failure is an explicit missing
  result; no molecule may be silently removed.

The calculation is a 5-ps-per-window PIMD2 response fingerprint, not a converged
full ARROW/PIMD8 free-energy calculation.

## Calibration locked in advance

Only the four calibration molecules may select Gaussian-process hyperparameters.
For each prior family, select from length scales {0.08, 0.12, 0.18, 0.25, 0.35}
and observation-noise multipliers {1, 2, 4} by minimum Gaussian negative log
predictive density. Each calibration molecule is held out in turn; its posterior
is conditioned only on the fixed 0.1, 0.5 and 0.9 observations and all twelve
initially hidden coordinates are scored. Ties use the smaller length scale and
then smaller noise multiplier. The residual amplitude is the root-mean-square
calibration residual, clipped to [1, 20] kcal mol^-1. The resulting numeric
settings and their input hashes must be committed in a dated calibration lock
before prospective sentinel responses are generated.

No architecture, kernel family, grid, prior interpolation, metric, budget or
schedule may change after calibration results are seen.

## Priors and posterior mechanics

- **Generic prior:** the calibration-set mean dense response curve.
- **SolvAI-conditioned prior:** the target molecule's already-frozen,
  molecule-held-out structure prediction at lambda 0.1, 0.5 and 0.9, extended to
  the fixed grid by PCHIP interpolation/extrapolation.
- **Generic covariance:** squared-exponential kernel with calibration-selected
  amplitude and length scale.
- **Molecule-conditioned covariance:** the same kernel family, multiplied by a
  bounded [0.75, 2.5] local scale computed deterministically from the absolute
  curvature of the structure prior.
- **Observation variance:** five-contiguous-block SEM squared, multiplied by the
  calibration-selected noise factor. Numerical jitter is `1e-10 * max(A^2, 1)`.
- **Integral:** the trapezoidal linear functional on the fifteen-point grid.

The inherited 0.1, 0.5 and 0.9 values are the common initial observations. Every
method receives identical observed values at a given total-window budget.

## Frozen policies and controls

Budgets are total observed windows {3, 4, 5, 7, 9, 12, 15}.

1. Prior only (generic and SolvAI-conditioned), before target observations.
2. Conventional PCHIP/trapezoidal integration of observed points.
3. Fixed population order:
   `0, 1, 0.3, 0.7, 0.05, 0.95, 0.2, 0.4, 0.6, 0.8, 0.75, 0.85`.
4. Uniform maximin refinement with lower-lambda tie breaking.
5. Twenty random schedules, seeds 91000--91019; all are reported and their mean
   is the random comparator.
6. Curvature-only schedule from the current structure-prior interpolant.
7. Generic Bayesian quadrature: generic mean/covariance, next point maximizing
   reduction in integral variance.
8. SolvAI-conditioned Bayesian quadrature: molecule-conditioned mean/covariance,
   the same variance-reduction acquisition and the same observation cost.
9. Oracle: after hypothetically revealing each candidate, chooses the point with
   smallest error to the dense integral. It is labelled non-deployable.
10. Full dense fifteen-window reference.

Stopping is the first budget at which the 90% posterior integral half-width is at
most 0.10 kcal mol^-1, with a minimum of 3 and maximum of 15 windows. Calibration
is scored at the actual stopping time as well as fixed budgets.

## Metrics and uncertainty

Primary target: absolute error of the reconstructed annihilation integral relative
to the dense fifteen-window trapezoidal integral. Secondary metrics: hidden-window
MAE/RMSE, maximum absolute local response error, signed integral bias, 90% integral
interval width and coverage, windows, 5-ps production time, bead-windows, nominal
bead-steps, measured wall time and GPU-hours. Simulation reconstruction and
experimental endpoint metrics remain separate.

Molecule-level paired intervals use 100,000 bootstrap resamples with seed
20260828. Random schedules are summarized across all twenty fixed seeds before a
molecule-level paired comparison.

## Interpretation rules

At total-window budgets 5 and 7, a positive result requires either:

1. SolvAI-conditioned Bayesian quadrature improves mean absolute integral error
   over the strongest deployable fixed, uniform, random, curvature-only and
   generic-BQ comparator by at least 0.10 kcal mol^-1, with the 90% paired
   bootstrap upper bound below zero and improvement on at least 6 of 8 sentinels;
   or
2. it reaches an integral MAE of at most 0.20 kcal mol^-1 with at least 20% fewer
   observed windows than every comparator, with 90% interval coverage between
   75% and 100% at the claimed stopping budget.

A neutral result has an absolute paired difference below 0.10 kcal mol^-1 with an
interval spanning zero. A negative result is a worsening of at least 0.10 kcal
mol^-1, severe undercoverage below 75%, or failure to beat a simpler schedule at
matched cost. No result at budget 9 or greater can rescue a failed primary
5/7-window comparison, although the complete frontier is reported.

If the prospective sentinel passes, Direction B survives and a separately frozen
multi-fidelity test may proceed. If it fails, Direction B is killed for this
protocol and Direction C is not launched. A favourable calibration molecule,
single sentinel, functional family, endpoint error or oracle result cannot rescue
failure.

## Already-known information and prohibited adaptation

Known before this freeze: all parent SolvAI results; Phase 1 aggregate and
molecule-level results; existing three-point responses and trajectory prefixes;
the fact that conditioned sparse interpolation sometimes beat a generic Gaussian;
the chosen molecules' identities and chemistry; and approximate historical
per-window runtimes. Unknown: all 144 new response values.

The experimental hydration labels are excluded from acquisition, reconstruction,
hyperparameter selection, stopping and molecule eligibility. No post-hoc molecule,
lambda, component or time-prefix selection is permitted. PIMD4/PIMD8 escalation,
new Hamiltonians and Tier-B remain outside this freeze.
