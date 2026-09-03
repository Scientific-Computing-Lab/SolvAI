# Hostile audit of the Active SolvAI no-go result

**Immutable scientific source:** `8fb984c2eb26d016c6b81cf488f88dc667ca9cd3`  
**Audit branch:** `active-solvai-v2-diagnostics`  
**Scope:** forensic verification only; no model tuning, new simulation, Tier-B access, or revision of the registered no-go decision.

## Executive finding

The central no-go result is reproducible from the raw Arbalest outputs. An independent implementation regenerated all 180 response means, standard deviations and five-block standard errors exactly; regenerated all 1,576 replay rows to floating-point precision; and reproduced every aggregate metric and paired bootstrap interval in the canonical Phase 2 files. The active policy therefore did not lose because of a lambda-index, sign, unit, cost-accounting, schedule-alignment, bootstrap-pairing or Gaussian-conditioning error.

The audit did, however, find an important interpretive problem in the non-deployable oracle. Its reported 0.337 and 0.068 kcal mol^-1 errors were obtained by selecting lambda locations and scoring them against the **same finite, noisy 5-ps curve**. Those values are optimistically selected diagnostics, not independent evidence that stable molecule-specific placement headroom exists. The cross-block analysis prospectively specified after this report will quantify how much, if any, survives held-out noise.

Three additional implementation/reporting defects were found. None can reverse the no-go decision:

1. The second success branch in `summarize_dense_sentinel.py` omitted its registered 75--100% coverage requirement. No active result reached the prerequisite 0.20-kcal mol^-1 MAE, so the corrected branch remains false.
2. The registered stopping rule was not emitted as a stopping-time result. Independent reconstruction across every sequential budget found **zero** molecule-policy paths with a 90% half-width at or below 0.10 kcal mol^-1; the smallest was 0.5474 kcal mol^-1. Thus the statement that no Bayesian method stopped is correct, but the original result file did not explicitly implement that analysis.
3. Nominal coverage was evaluated against a noisy dense integral using only posterior uncertainty, without dense-reference uncertainty or covariance from shared windows. The reported 0.750/0.875 values are descriptive coverage against the finite common pool, not validated frequentist calibration. This does not affect the point-error failure.

The registered scientific conclusion remains unchanged: the tested active method did not beat the strongest cost-matched comparator. The narrower claim that the in-sample oracle establishes molecule-specific placement headroom is suspended pending the independent-noise audit.

## 1. Raw response reconstruction

The independent audit script is `active_solvai/scripts/audit_no_go_pipeline.py`; its machine-readable output is `active_solvai/results/v2_diagnostics/stage1/stage1_audit.json`.

| Check | Result |
|---|---:|
| Dense response rows reconstructed | 180/180 |
| Stored frames per window | 51 in every case |
| Energy-file SHA-256 mismatches | 0 |
| Maximum response-mean discrepancy | 0 |
| Maximum response-SD discrepancy | 0 |
| Maximum five-block-SEM discrepancy | 0 |
| Maximum total-versus-component dH/dlambda discrepancy | 5.8e-9 kcal mol^-1 |

The component check sums every reported `dHdL_*` term and compares it with the raw `dHdL` column for every stored frame. The residual is only decimal-output round-off.

### Units, sign and thermodynamic direction

Arbalest's documented internal energy and native `.ene` response unit is kcal mol^-1 (`CommonLib/NaturalConstants.h`; `skills/arbalest/references/full-reference.md`, “Unit conventions”). Every audited XML transitions from `LIGSolvated` at lambda=0 to `Solvent` at lambda=1 and annihilates `LIG`. The raw integral is therefore the annihilation free energy. The analysis correctly reports hydration as

`Delta G_hyd = - integral_0^1 <dH/dlambda> dlambda`.

No extra kcal/kJ conversion is applied or required for native Arbalest output.

### Lambda mapping and quadrature

All 180 XML configurations contain the registered ordered grid

`0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0`.

For every row, `TIPoint` maps exactly to the recorded lambda; no mismatch was found. Independently reconstructed non-uniform trapezoidal weights sum to exactly 1.0 and reproduce every saved dense integral.

## 2. Reuse of the inherited windows

All 36 inherited observations (12 molecules times lambda 0.1, 0.5 and 0.9) resolve to the intended historical campaign directory, match the successful response inventory by molecule and lambda, use `TIPoint` 2, 6 and 12 respectively, and match the energy-file path recorded in the immutable manifest. No inherited window was duplicated, shifted or substituted.

The new configurations were copied from the corresponding molecule's inherited protocol and changed only in the intended case path, `TIPoint`, output path and disabled coordinate-trajectory output. The force field, molecule, solvent box, PIMD bead count, thermostat, barostat, time step and production length remain aligned.

## 3. Replay and canonical-number reproduction

The audit independently reimplemented interpolation, covariance construction, finite-grid Gaussian conditioning, integral variance, fixed/uniform/random/curvature/acquisition schedules, direct integration, oracle selection, molecule aggregation and paired bootstrapping.

| Artifact | Independent agreement |
|---|---:|
| Replay rows | 1,576/1,576 keys identical |
| Schedule mismatches | 0 |
| Largest prediction-field difference | 1.78e-14 |
| Largest aggregate-metric difference | 1.42e-14 |
| Largest paired-bootstrap difference | 1.11e-15 |

This exactly recovers the registered headline values:

| Windows | Active BQ MAE | Uniform-direct MAE | Paired change | 90% interval | Same-curve oracle MAE |
|---:|---:|---:|---:|---:|---:|
| 5 | 1.701248 | 1.152880 | +0.548369 | [-0.101072, +1.158993] | 0.336580 |
| 7 | 1.607515 | 1.091747 | +0.515768 | [-0.392137, +1.490642] | 0.068266 |

The final report, manuscript, Supplement and canonical JSON consistently round these values as documented.

## 4. Observation budgets and measured costs

Every policy begins with the same inherited three windows and receives exactly the registered total of 5 or 7 observed windows. Random schedules are first averaged across the 20 frozen seeds within each molecule and only then enter the molecule-paired comparison.

Window count, 5-ps production, bead-window count and nominal bead-step count are identical at a given budget. Because wall time varies slightly by molecule and lambda, measured wall time cannot be numerically identical across different schedules. The observed mean differences are negligible relative to simulation cost:

- At five windows, Active and uniform/generic use the identical lambda set and have identical mean measured wall time. Active differs from fixed by 0.241 s (0.029%).
- At seven windows, Active differs from uniform by -0.811 s (-0.087%), from generic by -0.571 s (-0.062%) and from fixed by +1.064 s (+0.115%).

Thus “identical cost” is exact in observed-window and simulated-time units and matched to within 0.12% in measured wall time. No result depends on a material cost imbalance.

## 5. Information-flow and leakage audit

The deployable schedules do not access dense prospective responses or experimental hydration labels:

- `active_solvai_bq` uses the frozen molecule-held-out structure response prior, a covariance calibrated on the four calibration molecules, and lambda-specific expected SEMs from those calibration molecules.
- `generic_bq` uses the calibration-set population mean/covariance and the same calibration-only expected-noise mechanism.
- `uniform_solvai_bq` is deterministic maximin refinement from the three common windows.
- `fixed_solvai_bq` uses the preregistered fixed order.
- `curvature_solvai_bq` uses only curvature of the frozen structure prior.
- Random schedules use the preregistered seeds and a stable molecule-ID hash.

Actual target-molecule SEMs enter posterior conditioning only **after** the corresponding window is observed. Candidate selection uses calibration-set expected SEMs, as required by amendment 001. Experimental endpoint labels never enter Phase 2 scheduling, conditioning, stopping or scoring.

The oracle is intentionally the sole exception: it reads the complete target curve and dense integral. It is correctly labelled non-deployable, but its original numeric headroom is not independently evaluated.

## 6. Was uniform sampling advantaged?

No hidden data or different observation budget was given to uniform sampling. At five windows, Active, uniform and generic BQ all select the same set `{0.1, 0.3, 0.5, 0.7, 0.9}`; `active_solvai_bq` and `uniform_solvai_bq` are therefore exactly identical. The much better `uniform_direct` result at that budget comes from its conventional PCHIP/trapezoidal estimator rather than a schedule advantage.

At seven windows, uniform and Active may select different locations, but both receive seven same-length PIMD2 windows and essentially identical measured cost. Reporting `uniform_direct` as the strongest comparator is scientifically fair: the freeze explicitly included direct interpolation, and the simpler estimator is allowed to win.

One display defect should be corrected in any later presentation: the Phase 2 frontier figure omits `uniform_direct`, although it is the strongest registered comparator at five and seven windows. The report table contains it, so no numerical conclusion is hidden, but the figure alone is incomplete. The original figure is preserved unchanged.

## 7. Gaussian conditioning, noise and intervals

The finite-grid Gaussian equations are implemented correctly. The posterior mean and covariance match an independent matrix implementation to machine precision. Observation variance is `noise_inflation * SEM^2`, consistent with treating the registered factor as a variance inflation factor. Numerical jitter is negligible and positive.

Important limitations:

1. Calibration NLPD scores noisy held-out 5-ps means against latent posterior variance without adding held-out measurement variance. This is a defensible latent-response objective only if the short-window mean is treated as the target; it is not evidence of calibrated prediction of an underlying converged response.
2. Original coverage compares a posterior interval for the latent integral with a dense target built from the same finite trajectory means. It omits dense-target uncertainty and shared-window covariance. The reported coverage is therefore descriptive, not a validated uncertainty-calibration claim.
3. The original replay did not save first-passage stopping rows. The independent audit evaluated every sequential budget from 3 through 15. The minimum 90% half-width was 0.547441 kcal mol^-1 and no path reached 0.10, so the no-stop claim is correct.

## 8. Bootstrap pairing and aggregation

Molecule pairing is correct. Candidate and comparator rows are indexed by identical molecule IDs before subtraction. The 20 random schedules are averaged within molecule before paired resampling, preventing random schedules from inflating the apparent sample size. The 100,000 resamples use the frozen seed 20260828. Independent recomputation recovers every 90% and 95% interval to at least 1.1e-15 kcal mol^-1.

Selecting the strongest comparator by observed point estimate is conservative against Active, although it introduces winner selection into the comparator label. It cannot create a false Active advantage and does not alter the negative result.

## 9. Concrete defects and disposition

| ID | Defect | Severity | Effect on registered no-go | Disposition |
|---|---|---|---|---|
| D1 | Condition 2 code omitted registered coverage gate | Logic defect | None: Active never reached the prerequisite 0.20 MAE | Documented; corrected evaluation is also false; original output preserved |
| D2 | Actual stopping-time rows were not emitted | Reporting/implementation omission | None: independent all-budget reconstruction finds zero stops | Documented; original output preserved |
| D3 | Coverage ignores noisy-reference uncertainty/shared-window covariance | Statistical interpretation defect | None on point-error failure | Coverage wording restricted to descriptive common-pool coverage |
| D4 | Main Phase 2 frontier omits winning `uniform_direct` curve | Figure defect | None: table and text report the comparator | Original preserved; diagnostic presentation will include it |
| D5 | Oracle selects and scores on the same noisy curve | Substantive diagnostic bias risk | Does not affect Active failure; invalidates an unqualified headroom inference | Quantified only after a separately committed cross-block protocol |

No registered file was overwritten. Original outputs and hashes remain intact under `active_solvai/results/phase2/`; audit reconstructions are isolated under `active_solvai/results/v2_diagnostics/stage1/`.

## Stage 1 conclusion

The Active SolvAI no-go is computationally and statistically secure with respect to its primary point-error comparison. The raw data, transformations, schedules, costs and paired intervals reproduce exactly. The only potentially material scientific correction concerns the oracle: its striking apparent advantage is in-sample and may reflect selection on trajectory noise. Until cross-block evaluation is complete, the evidence supports “the deployable policy failed” but not “stable molecule-specific quadrature headroom exists.”
