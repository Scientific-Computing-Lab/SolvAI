# V3 post-Gate-1 power and reference-reliability protocol

**Frozen before calculation:** 2026-09-03  
**Results already known:** the inherited v1/v2 results and the negative v3
Gate-1 identifiability result.  
**Purpose:** quantify—not authorize—the minimum independent design capable of
resolving the registered allocation effect under conservative inherited
heterogeneity and noise.

## Inputs

1. Molecule-level Active-SolvAI-BQ minus generic-BQ absolute integral-error
   differences on the eight v1 prospective dense molecules, at the registered
   five- and seven-window budgets.
2. The three unique complementary 2.5-ps dense-integral disagreements per
   molecule from the v2 independent-noise audit. Reverse-direction duplicates
   are not counted twice.
3. Measured v1 production throughput: 369.902 PIMD2 ps per GPU-hour.

These are development-exposed planning inputs. They do not establish v3 policy
performance.

## Effect and Monte Carlo design

The generic-BQ reference MAEs are 1.709526 kcal mol-1 at five windows and
1.114498 kcal mol-1 at seven windows. The minimum-useful alternative is a 20%
SolvAI reduction, giving true paired means of -0.341905 and -0.222900 kcal
mol-1. A 30% reduction is the fixed larger alternative; zero is the null.

For each budget, the eight observed paired differences are centered to zero.
Each Monte Carlo cohort resamples these molecule patterns with replacement and
adds the fixed alternative mean. Detection requires the upper endpoint of a
two-sided 95% paired t interval to be below zero. This is a tractable planning
surrogate for the final molecule-clustered bootstrap; it is not substituted for
the final analysis. Fifty thousand trials are run with base seed 20260903.
Type-I behavior is evaluated under the zero alternative.

Candidate molecule counts are 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128,
160, 192 and 256. Candidate production durations per independent stream-window
are 20, 50, 100, 200, 500, 1,000, 1,500, 2,000, 3,000, 4,000 and 5,000 ps.

## Independent-reference model

For each simulated cohort, molecule/partition absolute disagreements are
resampled from the 24 unique v2 values and scaled by `sqrt(2.5 / T)`, the
optimistic inverse-square-root trajectory-length assumption. The reliability
criterion is the master-freeze rule: the cohort 90th percentile of absolute
dense-integral disagreement must be no larger than half of the smaller
minimum-useful effect, 0.111450 kcal mol-1. The reported probability is the
fraction of Monte Carlo cohorts passing that rule.

Because the scaling assumption ignores persistent initialization and slow-mode
bias, any passing duration is a lower-bound planning estimate, not a validated
reference duration.

## Cost and decision

Production cost is `N * 15 windows * 2 streams * T / 369.902` GPU-hours.
Operational reservation is 1.11 times production cost. A grid point is
numerically adequate only if minimum-effect power is at least 0.80 at both
budgets, null type-I error is at most 0.075 at both budgets, and reference
reliability probability is at least 0.80.

Even a numerically adequate grid does not authorize simulation after the
negative Gate-1 result. A launch would additionally require a scientifically
credible reason that independent longer trajectories should reveal SolvAI
information absent from every frozen prefix, plus a separately committed
development-campaign freeze. If no affordable grid passes, v3 closes at the
quantitative gate and reports the smallest mathematical resolving design and
cost without launching it.
