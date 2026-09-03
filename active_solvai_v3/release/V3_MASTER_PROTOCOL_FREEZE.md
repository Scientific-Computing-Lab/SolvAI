# Active SolvAI v3 master protocol freeze

**Frozen:** 2026-09-03, before v3 difficulty-model fitting or held-out scoring  
**Branch:** `active-solvai-v3-effort-allocation`  
**Parent:** `770e13ce68ca80104557fd67a224abc5a2c44767`  
**Blueprint SHA-256:** `264f7f4afad7a9fa578992b33e6d5fefa83fe15f74be7fb9d85fd7e2f751b4de`

This document prospectively fixes the v3 causal question, information firewall,
zero-new-simulation analyses, launch gates and interpretation. It authorizes no
new molecular simulation. Later changes require a separately dated amendment
that states what results were already known.

## 1. Scientific question and estimand

The primary question is whether frozen molecule-conditioned SolvAI response
information reduces the physical effort needed to estimate a same-Hamiltonian
hydration free-energy integral when additional trajectory time is allocated
over a fixed 15-point alchemical grid.

The primary future causal contrast is:

> SolvAI-conditioned adaptive allocation minus a matched generic adaptive
> allocator that uses only lambda/protocol information and observations already
> revealed by the target simulation.

Generic adaptation versus equal time is a separate result. It cannot establish
an Active SolvAI contribution. The primary physical estimand is error and
uncertainty of the same-Hamiltonian dense integral at measured work. Experimental
hydration error is secondary and may be opened only after the primary simulation
analysis and interpretation are locked.

## 2. Fixed physical protocol

All primary policies retain lambda = {0.00, 0.05, 0.10, 0.20, 0.30, 0.40,
0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00}. They use the same
Hamiltonian, annihilation convention, PIMD2 settings, initialization rules,
quadrature, uncertainty estimator, mandatory initial floor, immutable chunk
size, hard per-window cap and total cost checkpoints. Policies may differ only
in the allocation of *additional* chunks. Lambda-location optimization is out
of scope.

For inherited Arbalest outputs, `SYSTEM dHdL` is read in native kcal mol-1,
lambda is sorted numerically, duplicate molecule-window records are forbidden,
the recorded annihilation convention is retained and the hydration sign is the
negative of the annihilation integral. Dense integration uses the fixed
trapezoidal weights of the registered grid. The t=0 production record is
excluded from time-series summaries. These conventions receive executable unit
tests before Gate-1 scoring.

## 3. Information firewall

At decision k, a deployable policy may access only fixed protocol metadata,
the target SMILES and deterministic structure representation, frozen SolvAI
response predictions, and chunks already revealed in its allocation stream.
It may not access future frames, unqueried chunks, the opposite reference
stream, a dense integral, experimental hydration labels, molecule-level future
failure status, prospective aggregate statistics or test-molecule results.

The independent reference stream is hidden until every action, stopping
decision and QC disposition for the allocation stream is immutable. Scoring is
performed A-to-B and B-to-A, then aggregated by molecule. Contiguous blocks of
one trajectory are never called independent replicas.

All exact identities and standardized equivalents remain in one fold. ARROW,
Tier-A, the historical PIMD2 probes and all twelve dense v1 molecules are
development-exposed. They cannot serve as the v3 prospective sentinel.

## 4. Gate-1 inherited development set

The full-grid Gate-1 set is fixed to the twelve inherited dense molecules:
Propane, Ethanol, Acetone, Pyridine, Cyclohexane, Octane, Octanol, AceticAcid,
Acetamide, DiMethylEther, EthylAcetate and Anthracene. Every molecule has the
same fifteen lambda values and one 5-ps PIMD2 stream per window. The 72 complete
three-point probe molecules may be used only for explicitly labelled secondary
diagnostics at lambda 0.1, 0.5 and 0.9; they cannot validate a full-grid policy.

The inherited stream is divided deterministically into ten consecutive 0.5-ps
blocks after dropping t=0. Nested diagnostic prefixes are 0.5, 1.0, 2.0 and
3.0 ps. The primary prefix is 1.0 ps; all four are reported as reliability
curves, with no favorable-length selection.

For every molecule-window, complementary later difficulty is computed only
from the non-overlapping final two 1-ps blocks (3.1--4.0 and 4.1--5.0 ps). The
primary stabilized target is

`log(0.5 * (mean_block_4 - mean_block_5)^2 + 0.25^2)`.

This is a realized complementary-block difficulty proxy, not an independent
variance truth. Stabilizers 0.10 and 0.50 kcal mol-1 are fixed sensitivity
analyses. Neither later block may enter a prefix feature.

## 5. Prefix diagnostics

At each frozen prefix the analysis computes, where estimable:

- response mean and sample variance;
- non-overlapping and overlapping batch-means variance rates;
- a Bartlett/Newey-West asymptotic variance rate with a deterministic
  data-length bandwidth;
- lag-1 autocorrelation, initial-positive-sequence integrated autocorrelation
  time and effective sample size;
- first-versus-second-half mean, variance and empirical-distribution distance;
- linear drift and rank drift;
- finite-estimator and unresolved-equilibration flags;
- available temperature/density diagnostics and measured/proportional cost.

Undefined short-prefix statistics remain explicit missing/flag values. They are
not silently clipped into favorable finite estimates. Reliability is reported
over prefixes and reasonable fixed block definitions.

## 6. Predeclared difficulty models

All scoring is leave-one-molecule-out. Every preprocessing transform and ridge
penalty is fit using training molecules only. The inner selection is grouped
leave-one-training-molecule-out over alpha in {0.1, 1, 10, 100, 1000}; ties
choose the larger penalty. Features are standardized using training rows only.
No tree or neural model is permitted at Gate 1.

The fixed model ladder is:

1. **Lambda/protocol:** lambda, lambda squared, endpoint indicator and fixed
   trapezoidal weight.
2. **Structure cold start:** lambda/protocol plus ten fixed RDKit descriptors:
   molecular weight, heavy-atom count, TPSA, MolLogP, H-bond donors, H-bond
   acceptors, rotatable bonds, ring count, fraction Csp3 and formal charge.
3. **SolvAI cold start:** lambda/protocol plus the exact fifteen frozen SolvAI
   response coordinates.
4. **Generic observed-only:** lambda/protocol plus the prefix diagnostics.
5. **Structure-conditioned observed:** generic observed-only plus the ten fixed
   descriptors.
6. **SolvAI-conditioned observed:** generic observed-only plus the fifteen
   frozen SolvAI response coordinates. This is the primary AI candidate at the
   identifiability gate.

The SolvAI model is compared primarily with generic observed-only; ordinary
structure is a diagnostic comparator, not a substitute causal baseline.
Response removal is model 4. The destructive alignment control jointly
permutes the complete fifteen-coordinate vector between training molecules,
using seeds 88031, 88032, 88033, 88034 and 88035. Test-molecule priors are never
permuted into training; the permutation is generated within each outer-training
set, and performance is summarized across all five fixed shuffles.

## 7. Gate-1 metrics and fixed interpretation

The primary prediction outcome is complementary log-difficulty at the 1-ps
prefix. Metrics are molecule-balanced MAE and RMSE, mean within-molecule
Spearman rank correlation over the fifteen windows and within-molecule
top-quartile difficulty AUROC. All molecule-level results are retained.

Paired uncertainty uses 100,000 molecule-clustered bootstrap resamples, seed
20260903. The primary difference is SolvAI-conditioned minus generic
observed-only molecule-mean absolute log error; negative favors SolvAI. The
aligned-minus-mean-shuffled and structure-conditioned-minus-generic contrasts
are also reported. Ninety- and 95-percent intervals are retained; the gate uses
the 90-percent interval because it is an early design-identifiability screen,
not a final claim.

SolvAI cold-start/conditioning is **identified as useful enough to inform a new
development design** only if all of the following hold at 1 ps:

1. SolvAI-conditioned MAE is at least 10% lower than generic observed-only and
   the 90% paired interval for the difference lies below zero;
2. aligned priors beat the mean shuffled control by at least 10% and the 90%
   paired interval lies below zero;
3. mean within-molecule Spearman correlation is at least 0.30; and
4. the direction does not reverse at the fixed 2-ps sensitivity prefix.

Generic diagnostic identifiability is recorded separately by the same 10%
and 90%-interval criterion versus lambda/protocol. Failure of the AI criterion
does not permit changing features, targets, prefixes, penalties or molecules.
A directional but uncertain result proceeds only to a prospective power
calculation. A null, reversed or shuffle-equivalent result closes the affected
cold-start claim unless a quantitatively powered independent development
design is justified from the frozen analysis.

No Gate-1 outcome is called prospective, independent-replica validation or a
deployable policy result.

## 8. Sequential replay contract

Before policy evaluation, the replay engine must:

- expose only the common initial chunk and then one immutable chosen chunk at
  a time;
- record the complete candidate state, utilities, tie-breaking seed, action,
  stop decision, model hash and cumulative cost at every step;
- reject attempts to read future chunks or the reference stream;
- reproduce byte-identical action ledgers from the same inputs and seed;
- pass synthetic tests in which the known hardest window and variance scaling
  are recovered without future-data access.

Inherited single-stream replay is implementation validation and development
diagnosis only. Scientific policy scoring requires genuinely independent
streams.

## 9. Future policy ladder

If and only if an independent development campaign passes its separate launch
gate, the fixed ladder will contain: equal-time allocation, random-chunk
allocation, naive maximum-SEM allocation, a competent generic shrinkage
allocator using target-simulation data only, a matched SolvAI-conditioned
shrinkage allocator, and a strictly development-only oracle. Every policy
receives identical initial chunks, eligible actions, caps, cost checkpoints,
quadrature and scoring references. A later `V3_POLICY_FREEZE.md` must specify
the one final model, shrinkage rule, uncertainty calculation and stopping rule
before a new sentinel is qualified.

## 10. Cost accounting

New physics cost, if separately authorized, includes equilibration,
initialization, production, failed jobs, retries, checkpoint/restart overhead
and policy orchestration. The ledger records simulated picoseconds,
force-evaluation-equivalent work where available, bead-windows, GPU seconds,
wall time, device and failure reason. A reserve must include at least 10% retry
contingency. Existing replay reports both equal nominal chunk cost and any
measured historical runtime; measured operational cost is primary for a future
sentinel.

## 11. Reference reliability, power and launch rules

An independent stream is not automatically a reference. Before any new
development or sentinel run, simulation-based power must include molecule
heterogeneity, two stream directions, noisy references, policy stochasticity
and the frozen molecule-clustered analysis. It must evaluate null, minimum-useful
and larger effects and report type-I behavior, expected interval width and
reference-gate probability.

The future reference gate is fixed in form: the 90th percentile of absolute
A-versus-B dense-integral disagreement must be no more than one half of the
minimum policy difference the design is powered to detect. A launch requires
at least 80% power for the minimum useful SolvAI increment, at least 80%
probability of satisfying the reference gate and an explicit storage/runtime
budget with 10% contingency.

The minimum useful prospective effect is fixed provisionally at either a 20%
reduction in independent integral MAE or normalized cost-error AUC versus the
generic allocator, or a 25% reduction in measured cost-to-target, with a paired
confidence interval excluding zero and no material loss of 80/90/95% coverage.
This threshold may be tightened, but not relaxed, in the development-campaign
freeze after realistic variance is quantified and before new trajectories.

**This master freeze authorizes zero new GPU simulation.** If Gate 1 cannot
calibrate a reliable design, the next artifact is a quantitative stop report.
If it supports a bounded campaign, a separate
`V3_DEVELOPMENT_CAMPAIGN_FREEZE.md` must name every molecule, identity hash,
seed, duration, chunk, retry rule, storage requirement, GPU ceiling and stop
condition and must be committed before launch.

## 12. Prospective interpretation

- Equal-time beaten by generic adaptation: evidence for generic adaptive
  effort allocation, not Active SolvAI.
- SolvAI-conditioned beaten by or indistinguishable from generic: no Active
  SolvAI contribution, regardless of performance versus equal time.
- Prospective SolvAI advantage satisfying the effect, interval, reliability and
  calibration gates: evidence sufficient to consider a predeclared scale-up.
- Unreliable references or inadequate power at the affordable ceiling: no run;
  report the exact minimum resolving design and cost.
- Experimental endpoint improvement cannot rescue a failed primary simulation
  result.

All failures, capped unresolved molecules and excluded runs remain in the
machine-readable ledgers. No favorable budget, family, molecule subset,
duration or target is selected after results are visible.
