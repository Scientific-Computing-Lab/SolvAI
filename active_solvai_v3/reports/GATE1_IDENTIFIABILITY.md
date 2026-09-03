# Gate 1: zero-new-simulation identifiability

## Decision

**The preregistered SolvAI identifiability gate did not pass.** Early target-
simulation diagnostics carry modest information about later complementary-block
difficulty, but the fifteen frozen SolvAI response coordinates add no stable
increment beyond those diagnostics. No new simulation is authorized by this
result.

This is a development-exposed measurement on complementary parts of one 5-ps
trajectory per molecule-window, not independent-replica validation.

## Frozen design

The analysis followed `release/V3_MASTER_PROTOCOL_FREEZE.md` and amendment 001,
committed at `4608c21` before model fitting. It used all 12 inherited dense
molecules, the fixed 15-point grid and four fixed prefix lengths. The primary
1-ps prefix was evaluated against a target built solely from the non-overlapping
3.1--4.0 and 4.1--5.0 ps blocks. Models were ridge regressions with nested
molecule-grouped penalty selection; every molecule and all 15 windows were held
out together. Five predeclared molecule-level response shuffles were included.

Command:

```bash
PYTHONPATH=active_solvai_v3/src uv run --project active_solvai_v3 \
  python active_solvai_v3/scripts/run_gate1_identifiability.py
```

## Primary results

| Model | Complementary log-difficulty MAE | Mean within-molecule rho | Top-quartile AUROC |
|---|---:|---:|---:|
| Lambda/protocol | 1.547 | 0.250 | 0.585 |
| Generic observed diagnostics | **1.427** | 0.248 | 0.573 |
| Structure-conditioned observed | 1.488 | 0.248 | 0.557 |
| SolvAI cold start | 1.676 | 0.147 | 0.589 |
| SolvAI-conditioned observed | 1.475 | 0.225 | 0.556 |

The generic diagnostic model improved MAE over lambda/protocol by 0.120 log
units (paired 90% CI -0.199 to -0.041), a 7.7% reduction. Its uncertainty
excluded zero, but it missed the frozen 10% materiality threshold; the generic
gate therefore did not formally pass.

The SolvAI-conditioned model was **worse** than the matched generic model by
0.0477 log units (90% CI -0.0171 to 0.1199), with only 6/12 molecules improved.
This is a -3.34% relative gain. The aligned model improved by only 0.00688
relative to the mean of the five shuffled controls; the paired difference was
-0.0102 (90% CI -0.0478 to 0.0310). Its mean within-molecule rank correlation
was 0.225, below the frozen 0.30 threshold.

The result is not a single-prefix accident:

| Prefix | Generic MAE | SolvAI-conditioned MAE | SolvAI minus generic |
|---:|---:|---:|---:|
| 0.5 ps | 1.438 | 1.466 | +0.028 |
| 1.0 ps | 1.427 | 1.475 | +0.048 |
| 2.0 ps | 1.430 | 1.438 | +0.008 |
| 3.0 ps | 1.444 | 1.460 | +0.016 |

Positive differences favor the generic model. Results were also directionally
unchanged under the 0.10 and 0.50 kcal mol-1 target stabilizers: the SolvAI
increment remained +0.062 and +0.025 log units at the primary prefix.

## What was identifiable

At 1 ps, the strongest raw relationships with complementary difficulty were
the five-batch variance rate (Spearman rho=0.362), naive SEM or sample variance
(rho=0.337), Newey--West variance rate (rho=0.326) and overlapping-batch rate
(rho=0.320). These modest correlations support the physical premise of generic
adaptive sampling, but not a molecule-conditioned SolvAI increment. At 3 ps,
generic rank correlation rose to 0.354 and top-quartile AUROC to 0.699, while
the SolvAI-conditioned MAE remained worse than generic.

## Interpretation

The response coordinates that improve the experimental SolvAI endpoint do not,
on these data, identify which fixed-grid windows will have large later sampling
disagreement. Cold-start SolvAI features alone were worse than the lambda-only
baseline by 0.129 log units (90% CI 0.051 to 0.225). Ordinary structure also did
not improve the generic observed model. This pattern argues against spending
GPU time merely to fit the frozen v3 candidate on longer versions of the same
twelve molecules.

The inherited data remain too short and non-independent to prove that a useful
effect is absent in the wider chemical population. The next step is therefore
not a simulation launch but the separately frozen power/reference-reliability
calculation: it must quantify the minimum independent design that could resolve
the registered 20% SolvAI-over-generic effect and its physical-reference cost.

## Machine-readable evidence

- `results/gate1/gate1_prefix_diagnostics.parquet`
- `results/gate1/gate1_oof_predictions.parquet`
- `results/gate1/gate1_model_metrics.csv`
- `results/gate1/gate1_paired_comparisons.csv`
- `results/gate1/gate1_diagnostic_reliability.csv`
- `results/gate1/gate1_canonical_metrics.json`
- `results/gate1/gate1_artifact_manifest.json`
