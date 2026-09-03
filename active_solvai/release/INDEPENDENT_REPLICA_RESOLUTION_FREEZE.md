# Independent-replica resolution experiment: prospective freeze

**Frozen before the formal power calculation and before any new simulation:** 2026-09-03 UTC  
**Immutable Active SolvAI no-go:** `active-solvai` at `8fb984c2eb26d016c6b81cf488f88dc667ca9cd3`  
**Immutable independent-noise diagnostic:** `active-solvai-v2-diagnostics` at `f2cc60fb73416e0e417c97c600d5abd00868031c`  
**Execution branch:** `active-solvai-v2-independent-replicas`

## Scientific status and sole question

The failed Active SolvAI acquisition policy and its registered conclusions remain immutable. This is neither a rescue of that policy nor a new experimental-endpoint test. The sole question is whether *stable molecule-specific lambda-placement headroom* exists when finite-trajectory response noise is reduced with independently initialized, longer PIMD2 replicas.

No model, prior, covariance function, acquisition rule, lambda grid, estimator, budget or interpretation threshold may be tuned using these replicas. No Tier-B data, experimental hydration endpoint or new PIMD8 calculation may be used.

## Inputs fixed before power calculation

The power calculation may read only the already development-exposed diagnostic artifacts below.

| Input | SHA-256 |
|---|---|
| `active_solvai/results/v2_diagnostics/oracle_independent_noise/cross_block_predictions.parquet` | `e6492f5754aab9d4724e10c9b5e2d08863e669beb9b252a4f68a3e0b88a715fa` |
| `active_solvai/results/v2_diagnostics/oracle_independent_noise/molecule_level_metrics.csv` | `7be6b7c3580b8386e6f31395463a9a6d6723d41c1590f2d714a066dd4fff3679` |
| `active_solvai/results/v2_diagnostics/oracle_independent_noise/dense_integral_reproducibility_by_split.csv` | `feec7765a4705ca9390e87fd83eeb29847539cbc18f2aaa16f08bc7c71c3df0e` |
| `active_solvai/results/v2_diagnostics/oracle_independent_noise/oracle_schedule_stability.csv` | `09a30eb8afacd2220544e4f6f07163e8655fa2b6975e3a2501669d80eda2cdeb` |
| `active_solvai/results/v2_diagnostics/oracle_independent_noise/canonical_metrics.json` | `c0f01f65f3cfae4230d01fb536b7f375a1a0db29742a5023870ed53022cec2f4` |
| `active_solvai/results/phase2/dense_responses_prospective.parquet` | `4ce191a2ebfabb10dfcf5e11f98fef91b0bd5b2993f066ab47cfa9cc261c1097` |
| `active_solvai/results/phase2/dense_replay_predictions.parquet` | `a0dbc9b61c7f54848be8e56dcc5c22d12d511839c90cdbfe44900574112bf1fd` |
| `active_solvai/configs/dense_sentinel_v1.json` | `fc6afd2a4f8a255a641228df624769e6ed437dc4cc69dc244480a33c2936594d` |
| `active_solvai/simulations/dense_pimd2/manifest.csv` | `9c46ad5b78af027f3c82d84ff300c1ea839fbd13650eba704bbfef95f5871358` |

The script must assert all hashes before calculation. The diagnostic's observed 2.5-ps complementary-half dense-integral MAE (1.86575 kcal mol-1), oracle-minus-uniform paired effects and reversal schedules are known inputs, not results of this experiment.

## Pre-specified power and precision gate

### Estimand

At each total-window budget `b` in `{5, 7}`, define the paired molecule effect

`d_i,b = |G_oracle,i,b - G_dense,i| - |G_uniform-direct,i,b - G_dense,i|`.

Negative values favor oracle placement. The meaningful alternative is a population mean of `-0.30 kcal mol-1`; smaller differences are not used to justify simulation.

### Noise and heterogeneity model

1. Use the six registered cross-block directions per molecule and budget. Average each exact reversal pair, yielding three partition effects per molecule.
2. The eight molecule means define the empirical between-molecule pattern after centering it and shifting its mean to exactly `-0.30 kcal mol-1`. This preserves observed heterogeneity without treating the noisy observed positive mean as truth.
3. Within-molecule residuals are the three partition effects minus their molecule mean. Their pooled variance is the empirical 2.5-ps effective-half noise. For a proposed independent-replica length `T`, scale this variance by `2.5/T`, the observed inverse-time law. This is explicitly an idealized precision projection; departures from inverse-time behavior can only make the campaign less reliable.
4. Simulate two independent transfer directions per molecule and average them exactly as the planned analysis will. For the fixed eight-sentinel design, retain all eight centered molecule effects. For alternative sample sizes, sample the centered eight-molecule effects with replacement, independently of the simulated response noise.
5. Use 100,000 Monte Carlo experiments per design and seed `20260904`.

### Detection event and power threshold

A simulated experiment detects material headroom at a budget only if all three conditions hold:

- mean `d <= -0.20 kcal mol-1`;
- the one-sided 90% paired Student interval has an upper bound below zero; and
- at least 75% of molecules have `d < 0`.

The Student interval is used only for prospective power because evaluating 100,000 nested 100,000-resample bootstraps is impractical. The actual result will use the registered molecule-clustered bootstrap.

The proposed eight-molecule, 50-ps design is adequately powered only if detection probability is at least 0.80 at **both** five and seven windows.

### Dense-reference precision

For each of the three complementary partitions and molecule, the observed absolute 2.5-ps half-integral difference is treated as the magnitude of a zero-mean normal difference at that molecule's observed scale. Independent 50-ps replica differences are simulated by multiplying by `sqrt(2.5/T)` and a half-normal draw normalized to unit mean. A simulated dense reference passes when its across-molecule MAE is no more than 0.50 kcal mol-1 and its median absolute difference is no more than 0.30 kcal mol-1. The design requires at least 0.80 probability of passing both criteria.

An additional conservative sensitivity uses the unscaled empirical between-molecule variance and doubles the within-molecule variance; it is reported but is not substituted after results are seen.

### Launch rule and alternative design search

Launch the registered eight-molecule campaign only if the primary model reaches 0.80 power at both budgets and 0.80 dense-reference reliability. Otherwise **no simulation is launched**.

If it fails, search the fixed grid:

- molecules `n in {8, 12, 16, 24, 32, 48, 64, 96}`;
- independent production per replica `T in {50, 75, 100, 150, 200} ps`;
- two replicas, all 15 windows.

The smallest adequately powered alternative is the passing design with the least projected production GPU-hours; ties use fewer molecules, then shorter trajectories. This is a design recommendation only. Molecules beyond the original eight would require a separate chemistry-first prospective selection freeze and are not authorized here.

Measured throughput is frozen from the original prospective pool: 480 ps in 1.2976373865 GPU-hours (`369.902 ps GPU-hour-1`). Projected production cost is `n * 15 * 2 * T / 369.902`; the operational reservation is 1.11 times this value. No alternative simulation is launched.

The power script, full simulation draws summary, passing-design grid and report are saved under `active_solvai/results/v2_independent_replicas/power/` and `active_solvai/reports/`.

## Authorized campaign, conditional on passing the gate

### Molecules and lambda grid

Exactly these eight sentinels, with no substitution or exclusion:

`Cyclohexane, Octane, Octanol, AceticAcid, Acetamide, DiMethylEther, EthylAcetate, Anthracene`.

Every molecule uses all 15 registered lambdas:

`0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0`.

### Frozen simulation protocol

- Force field: ARROW 933; solvent: water; 24-Angstrom box.
- Temperature: 298 K; pressure: 1 bar; PIMD beads: 2.
- Time step: 0.002 ps; fast-force divider: 16.
- Frozen minimization: 250 steps.
- Frozen equilibration: 500 steps (1 ps), PILE thermostat, matching inherited protocol.
- Production: exactly 25,000 steps = 50 ps per replica, PILE thermostat and Berendsen barostat.
- Energy output: every 50 steps (0.1 ps; 501 stored values including the initial production record when emitted by the engine).
- Coordinate trajectory output remains disabled unless the engine requires a minimal restart checkpoint; energy and checkpoint outputs are retained.
- Two independent replicas per molecule and lambda; total 240 trajectories, 12 ns PIMD2 production and 480 bead-windows.

The starting molecular/solvent coordinates and all non-random physical settings are inherited unchanged from the corresponding registered dense-sentinel configuration. Independence is introduced by both velocity and PILE thermostat seeds.

### Deterministic seed schedule

Molecules use their order above (`m = 0..7`), lambdas their order above (`l = 0..14`) and replicas `r = 0,1`.

- velocity seed: `51000001 + 100000*r + 1000*m + l`;
- PILE thermostat seed: `61000001 + 100000*r + 1000*m + l`.

The seed is inserted explicitly in every initialization and PILE thermostat block. A configuration validator must prove all 480 stochastic seeds are distinct and match this formula. A retry resumes the same configuration/checkpoint; it never silently changes a seed. If a fresh restart is scientifically required, it is recorded as a failed run and uses predeclared offset `+1000000 * retry_index` for both seeds. The original attempt remains in the ledger.

### QC and result embargo

Before production, one short technical smoke test may validate parsing, output cadence, temperature, finite energies, restartability and GPU use; it may not calculate or display oracle/comparator performance.

All 240 registered trajectories must be accounted for. Automated QC before unblinding checks only:

- successful engine exit or documented retry/failure;
- exact configuration and input hashes;
- correct molecule, lambda, bead count, step count and seeds;
- finite energies and dH/dlambda values;
- expected output cadence and no unexplained missing frames;
- absence of catastrophic temperature/energy instability;
- complete compute/failure ledger.

No molecule or window is silently excluded. Oracle, Active, comparator, dense-integral or schedule-stability performance is not computed or inspected until both replicas pass all prespecified QC for every trajectory. If irrecoverable failures remain, all intended rows are reported and the registered performance analysis is not represented as complete.

## Frozen analysis

### Independent-replica transfer

At 50 ps, perform both directions separately:

1. Replica A selects oracle lambda locations; replica B supplies the selected response values and the dense 15-window target.
2. Replica B selects; replica A evaluates.

Only after reporting each direction are they averaged within molecule. The oracle uses the unchanged greedy oracle BQ selector. Active SolvAI BQ, generic BQ, uniform direct, fixed direct and the 20 registered random schedules use their original frozen definitions. Every method receives exactly five or seven 50-ps windows in the evaluation replica. The three inherited lambda locations remain `{0.1, 0.5, 0.9}` conceptually, but all values now come from the corresponding new independent replica; no old 5-ps response is mixed into this result.

No experimental hydration endpoint is read.

### Reported quantities

For every molecule, direction, method and budget save selected lambda indices/values, dense target, estimate, signed error, absolute error, posterior uncertainty and interval inclusion where defined.

Report:

- dense-integral A-versus-B signed differences, MAE, RMSE, median absolute difference, correlation and molecule-level values;
- each transfer direction and their within-molecule average;
- oracle-set Jaccard similarity between A and B, excluding the common inherited three;
- MAE, RMSE and molecule-level errors for every method at five and seven windows;
- oracle-minus-comparator paired molecule differences and fraction improved;
- 90% and 95% percentile intervals from 100,000 molecule-clustered bootstrap resamples, seed `20260904`;
- random schedules averaged across the 20 frozen seeds within molecule before aggregation;
- descriptive fixed-prefix analyses at 5, 10, 20, 30 and 50 ps, always using the prefix from each independent production replica.

The 50-ps result is primary. Prefixes may diagnose convergence but cannot replace, select or redefine it.

### Frozen reliability criterion

The dense reference is reliable only if the eight-molecule A-versus-B dense-integral MAE is at most 0.50 kcal mol-1 **and** its median absolute difference is at most 0.30 kcal mol-1. Failure of either threshold invokes decision C unless decision A is established by a clear absence/instability of placement advantage despite the noisy reference.

### Frozen decision rule

Let `delta_b` be the mean paired `oracle - uniform-direct` absolute-error difference after averaging the two directions within each molecule. Negative values favor oracle placement. Let `U90_b` be its molecule-clustered 90% upper confidence limit, `q_b` the fraction of molecules improved and `J_b` the mean A/B Jaccard similarity of added locations.

Apply the following hierarchy and return exactly one decision:

#### B. STABLE PLACEMENT HEADROOM ESTABLISHED

Return B only if the dense-reference reliability criterion passes and, at least one of five or seven windows:

- `delta_b <= -0.20 kcal mol-1`;
- `U90_b < 0`;
- `q_b >= 0.75`; and
- `J_b >= 0.40`;

while the other budget is not materially contradictory (`delta > +0.10` with a 90% lower confidence bound above zero). This authorizes designing, not running, an Active SolvAI v2 program to learn placement rules from a larger dense-curve corpus. It does not authorize Tier-B, PIMD8 escalation or an experimental-endpoint claim.

#### A. HEADROOM REFUTED

If B fails and the dense reference passes, return A when either:

- both budgets have `delta_b >= -0.10 kcal mol-1` and 90% lower confidence bounds above `-0.30 kcal mol-1`; or
- the placement sets remain unstable (`J_b < 0.40` at both budgets), even if a noisy point estimate appears favorable.

Close this Active SolvAI formulation. Do not launch further PIMD, Tier-B or acquisition tuning.

#### C. STILL INCONCLUSIVE AFTER 50 PS

Return C when neither A nor B is established because the dense reference fails its reliability criterion or molecule-level heterogeneity/uncertainty prevents either conclusion. C is a terminal practical no-go for this PIMD2 route; no trajectory extension is requested.

## Records, manuscript and immutability

All configurations, raw outputs, retries, hashes, wall/GPU time, failures, force-evaluation equivalents, bead-windows, predictions, metrics, plots and decisions are stored only under the new v2-independent-replica namespaces. The immutable prior branches and artifacts are never overwritten.

If the power gate fails, the experiment terminates before simulation. The required output is then the frozen power analysis, the smallest adequately powered alternative design and exact projected GPU cost; no A/B/C experimental decision is fabricated.

If the campaign runs, update the Active SolvAI report/manuscript/Supplement only with clearly labelled resolution-experiment notes. This experiment cannot revise the original policy no-go. Full tests, manifests, security checks, PDF inspection, `git diff --check`, a clean commit and branch push are required before handoff.
