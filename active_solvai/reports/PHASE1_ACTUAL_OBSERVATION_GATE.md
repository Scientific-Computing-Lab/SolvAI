# Phase 1 actual-observation gate

## Decision

**The preregistered experimental-endpoint gate is negative.** On the 72 molecules with complete 5 ps PIMD2 observations at λ=0.1, 0.5 and 0.9, the actual-minus-predicted response correction increased repeated-partition MAE from 0.186374 to 0.189856 kcal mol⁻¹. The molecule-paired candidate-minus-baseline change was +0.003482 kcal mol⁻¹ (95% bootstrap CI +0.000960 to +0.006093). Only 31.9% of molecules improved. The candidate was also indistinguishable from the mean shuffled-residual control (difference +0.000580; 95% CI -0.001191 to +0.002526).

The sign was unfavourable in four of five repeated partitions and effectively null in the fifth. The fixed parent partition was also worse (0.180392 versus 0.182904). This satisfies the frozen negative endpoint criterion and cannot be rescued by a favourable post-hoc λ subset or chemistry.

## Matched fixed-partition results

| method | n | mae | rmse | median_absolute_error |
| --- | --- | --- | --- | --- |
| P1-A_frozen_solvai | 72 | 0.1804 | 0.2874 | 0.1075 |
| P1-B_predicted_response | 72 | 0.1831 | 0.2894 | 0.1204 |
| P1-C_actual_response | 72 | 0.1834 | 0.2919 | 0.1096 |
| P1-D_actual_minus_predicted | 72 | 0.1829 | 0.2909 | 0.1121 |
| P1-H_mean_shuffled_residual | 72 | 0.1823 | 0.2895 | 0.1136 |

## Five repeated partitions

| method | mean_mae | sd_mae | min_mae | max_mae |
| --- | --- | --- | --- | --- |
| P1-A_frozen_solvai | 0.1864 | 0.0050 | 0.1806 | 0.1935 |
| P1-H_mean_shuffled_residual | 0.1893 | 0.0040 | 0.1839 | 0.1928 |
| P1-D_actual_minus_predicted | 0.1899 | 0.0040 | 0.1836 | 0.1934 |
| P1-C_actual_response | 0.1907 | 0.0042 | 0.1841 | 0.1939 |
| P1-B_predicted_response | 0.1916 | 0.0057 | 0.1849 | 0.1971 |

## Destructive control

All three observed response values were permuted jointly across molecules within each outer training and test fold for five preregistered seeds. The aligned residual was not better than the mean shuffled control. This falsifies the claim that the present 5 ps PIMD2 residual adds stable molecule-specific endpoint information beyond frozen SolvAI.

## Observation duration

The primary residual correction was repeated at sequential 0.5, 1.0, 2.0, 3.5 and 5.0 ps prefixes. Longer prefixes reduced variability but did not produce an endpoint gain over frozen SolvAI. No future frames were selected.

## Sparse response reconstruction

SolvAI-conditioned Gaussian interpolation often reduced hidden-point error relative to a generic population Gaussian. The largest post-result descriptive gain among the predeclared two-point subsets was `0p5_0p9`: 10.527 to 4.679 kcal mol⁻¹, paired difference -5.848 (95% CI -7.100 to -4.598). However, hidden-point errors remained 2.4 kcal mol⁻¹ or larger and nominal 95% intervals were very wide. These three-point held-point diagnostics are not dense-curve reconstruction. Under the prospective freeze they justify only a bounded dense same-Hamiltonian sentinel test, not a reconstruction claim or endpoint rescue.

## Numerical integration

Direct integration of the three actual short-window observations followed by fold-local affine calibration produced a five-repeat mean experimental MAE of 1.426 kcal mol⁻¹. The best posterior-integral variant remained 0.978 kcal mol⁻¹. The inherited three λ points are therefore not an adequate quadrature rule for this endpoint.

## Direction decisions

1. **Direction A — empirical residual correction:** failed the prospectively frozen endpoint criterion.
2. **Direction B — molecule-conditioned Bayesian quadrature:** hidden-point interpolation gains passed the limited diagnostic criterion for several predeclared subsets, but no compatible dense same-Hamiltonian population exists locally, intervals are broad, and the integral is inaccurate. This is a conditional go to one bounded, prospectively frozen dense sentinel acquisition; it is not a reconstruction claim.
3. **Direction C — adaptive multi-fidelity allocation:** held until that dense sentinel establishes whether a useful reconstruction model exists.

No new MD/PIMD calculation was used in this gate. This is a scientific null for endpoint correction with the inherited 5 ps, PIMD2, three-window protocol—not evidence that all possible active solvation calculations lack value. A separately frozen dense sentinel can test the surviving reconstruction hypothesis, but it cannot alter the failed endpoint decision.

## Reproducibility

- Freeze commit: `a0dd986`
- Command: `active_solvai/.venv/bin/python active_solvai/scripts/run_phase1_gate.py`
- Summary command: `active_solvai/.venv/bin/python active_solvai/scripts/summarize_phase1.py`
- Machine-readable endpoint, response, reconstruction and integration tables are in `active_solvai/results/phase1/`.
