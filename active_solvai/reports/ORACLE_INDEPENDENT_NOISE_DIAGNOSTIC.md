# Independent-noise audit of oracle placement headroom

**Pre-specified post-hoc conclusion: C. INCONCLUSIVE DUE TO RESPONSE NOISE.**

This diagnostic was frozen at commit `d0167dcdfc4298a3b0a0ffbfa0be1b87bc8f2be3`, before cross-block performance was calculated. It reuses the eight development-exposed prospective sentinels and their existing 5-ps PIMD2 trajectories. It is not new prospective molecular validation, does not reopen the failed Active SolvAI policy, and uses no new simulation.

## What was tested

Each 51-frame lambda trajectory was divided into four fixed contiguous blocks (13, 13, 13 and 12 frames). All six ordered two-block selection/two-block evaluation assignments were used. For the non-deployable oracle, lambda locations were selected only from the selection blocks. The selected locations were then populated with response values from the non-overlapping evaluation blocks and scored against a 15-window dense integral built only from those evaluation blocks. Reversing the assignments tests whether the chosen locations survive independent trajectory noise.

Active, uniform, fixed and generic-BQ policies were evaluated against the identical held-out-block dense targets. Uniform-direct and fixed-direct estimators were retained as conventional secondary controls. All methods received the same number of evaluation-side observations at each budget. Performance was averaged over all six splits within each molecule before molecule-level aggregation and bootstrap resampling.

## Dense-integral reproducibility

The short trajectories do not provide a stable dense target when divided into non-overlapping halves.

| Quantity | Result (kcal mol^-1) |
|---|---:|
| Complementary-half dense-integral MAE | 1.866 |
| 90% molecule-clustered interval | [1.418, 2.312] |
| 95% molecule-clustered interval | [1.329, 2.385] |
| Median absolute half difference | 1.382 |
| RMSE | 2.276 |

The prespecified reliability gate was 0.50 kcal mol^-1 mean or 0.30 kcal mol^-1 median. Both were exceeded by a wide margin. These are non-overlapping blocks from the same trajectory, not independently initialized replicas, so even this large disagreement can understate some slow-mode uncertainty.

Molecule-level mean complementary-half differences ranged from 0.832 kcal mol^-1 for cyclohexane to 3.011 kcal mol^-1 for octanol. The instability is therefore not confined to one chemical outlier.

## Cross-fitted replay

| Method | 5-window MAE | 90% interval | 7-window MAE | 90% interval |
|---|---:|---:|---:|---:|
| Cross-fitted oracle BQ | 1.732 | [1.270, 2.224] | 1.681 | [1.174, 2.199] |
| Active SolvAI BQ | 1.973 | [1.519, 2.428] | 2.122 | [1.544, 2.753] |
| Uniform SolvAI BQ | 1.973 | [1.519, 2.428] | 1.774 | [1.449, 2.114] |
| Fixed SolvAI BQ | 3.073 | [2.205, 4.066] | 2.559 | [1.773, 3.393] |
| Generic BQ | 1.800 | [1.298, 2.363] | 1.273 | [0.920, 1.627] |
| Uniform direct | 1.301 | [1.015, 1.617] | 1.272 | [1.000, 1.547] |
| Fixed direct | 1.835 | [1.479, 2.196] | 1.764 | [1.260, 2.277] |

Against uniform direct, the cross-fitted oracle is worse by +0.431 kcal mol^-1 at five windows (90% interval -0.281 to +1.138) and +0.409 kcal mol^-1 at seven windows (-0.233 to +1.075). It improves only three of eight molecule-level split averages at either budget. The intervals include both a material benefit and material harm, so the diagnostic cannot support either stable headroom or a precise null.

## How much oracle headroom survived?

| Windows | Original same-curve headroom | Cross-fitted headroom | Survival ratio |
|---:|---:|---:|---:|
| 5 | +0.816 | -0.431 | -0.528 |
| 7 | +1.023 | -0.409 | -0.399 |

Headroom is defined as `uniform-direct MAE - oracle MAE`; positive values favor the oracle. Thus none of the original point-estimate headroom survives held-out-block scoring. The sign reverses. This directly confirms that the original 0.337/0.068 values were strongly optimistically biased by using the same noisy curve for subset selection and evaluation.

That reversal alone does not meet the stricter predeclared criterion for conclusion A, because the molecule-clustered intervals remain too broad to rule out a 0.10-kcal mol^-1 benefit. The appropriate conclusion is therefore C rather than a definitive assertion that all placement headroom is absent.

## Lambda-selection stability

The selected added lambda locations are unstable under reversal of the trajectory blocks:

| Budget | Mean reversal-pair Jaccard | Mean all-pairs Jaccard |
|---:|---:|---:|
| 5 windows (2 added) | 0.139 | 0.189 |
| 7 windows (4 added) | 0.238 | 0.282 |

The same molecule therefore usually receives a substantially different oracle subset when the non-overlapping half used for selection changes. No single lambda dominates sufficiently to rescue the molecule-specific placement interpretation; selection frequencies are distributed across the 12 eligible added positions. Complete schedules and frequencies are available in `crossfit_oracle_schedules.csv`.

## Molecule-level result

| Molecule | Oracle 5 | Uniform 5 | Oracle 7 | Uniform 7 |
|---|---:|---:|---:|---:|
| Acetamide | 3.190 | 0.859 | 2.779 | 0.648 |
| Acetic acid | 0.792 | 1.525 | 0.946 | 1.275 |
| Anthracene | 1.939 | 0.960 | 1.423 | 0.931 |
| Cyclohexane | 1.255 | 0.746 | 0.889 | 0.792 |
| Dimethyl ether | 0.949 | 1.830 | 0.689 | 1.890 |
| Ethyl acetate | 0.975 | 2.312 | 1.080 | 1.818 |
| Octane | 2.245 | 1.266 | 2.932 | 1.018 |
| Octanol | 2.512 | 0.912 | 2.712 | 1.807 |

Values are mean absolute errors across all six fixed block assignments, in kcal mol^-1. The mixed molecule-level signs and broad intervals reinforce the noise diagnosis; no family or molecule was selected post hoc to alter the conclusion.

## Exact conclusion

### C. INCONCLUSIVE DUE TO RESPONSE NOISE

The existing 5-ps trajectories cannot support a reliable independent evaluation of molecule-specific oracle placement. Same-curve oracle selection is severely optimistic: its apparent advantage reverses on held-out blocks, selected subsets are unstable, and independent half-curve dense integrals disagree by substantially more than the target effect. Yet with only eight molecules and 2.5-ps evaluation halves, uncertainty is too broad to rule out a material true placement advantage.

This finding does **not** change the immutable Active SolvAI no-go. The deployable policy still failed its original prospective comparison. It changes only the interpretation of the oracle: the original replay did not establish stable headroom.

## Smallest resolving experiment—not authorized or launched

A credible resolution requires independent trajectory replicas, not further slicing of the same 5-ps trace.

- **Molecules:** the same frozen eight sentinels: Cyclohexane, Octane, Octanol, AceticAcid, Acetamide, DiMethylEther, EthylAcetate and Anthracene.
- **Windows:** all 15 registered lambda values for every molecule. A dense common pool is necessary because oracle selection cannot be assessed against unobserved locations.
- **Replicas:** two independently initialized replicas per molecule-window, with independent velocities and thermostat/random seeds.
- **Length:** 50 ps production per replica after the existing protocol's equilibration and QC, giving 100 ps total production per molecule-window.
- **Scale:** 240 trajectories, 12 ns aggregate PIMD2 production and 480 bead-windows.
- **Measured-cost projection:** the original prospective campaign delivered 480 ps in 1.2976 GPU-h. Linear scaling gives approximately 32.4 GPU-h of production-equivalent work; reserve **36 GPU-h and 40 h wall time on one RTX 3090** for equilibration, I/O and contingency.
- **Frozen analysis:** replica A selects and replica B evaluates, then reverse; use the same 15-point grid, priors, policies, budgets and molecule-clustered comparisons. No adaptive extension or favorable molecule selection.

The 50-ps choice is grounded in the observed 1.866-kcal mol^-1 half-curve disagreement: under ideal inverse-square-root scaling, increasing an effective 2.5-ps half to 50 ps would reduce it to approximately 0.42 kcal mol^-1, below the prespecified 0.50 reliability gate. This is a resolution experiment, not a justified scale-up, and it must not be launched without approval.

## Reproducibility

Run:

```bash
uv run --project active_solvai python active_solvai/scripts/run_oracle_independent_noise_diagnostic.py
```

Machine-readable outputs are under `active_solvai/results/v2_diagnostics/oracle_independent_noise/`; diagnostic figures are under `active_solvai/figures/v2_diagnostics/`. Input hashes and the exact conclusion rule are in the committed freeze.
