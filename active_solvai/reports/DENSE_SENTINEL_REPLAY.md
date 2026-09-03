# Prospective dense PIMD2 sentinel replay

**Frozen decision: NEGATIVE.**

This is a prospective test of short-window, same-Hamiltonian PIMD2 response reconstruction. It is not a test of experimental endpoint improvement, full PIMD8 convergence or blind chemistry.

## Primary frozen comparisons

| Total windows | Active MAE | 90% coverage | Strongest comparator | Comparator MAE | Paired difference | 90% CI | Improved | Pass |
|---:|---:|---:|---|---:|---:|---:|---:|:---:|
| 5 | 1.701 | 0.750 | uniform_direct | 1.153 | +0.548 | [-0.101, +1.159] | 0.250 | no |
| 7 | 1.608 | 0.875 | uniform_direct | 1.092 | +0.516 | [-0.392, +1.491] | 0.375 | no |

## Accuracy–cost frontier

| Method | Windows | Integral MAE | Hidden-curve MAE | 90% coverage | Mean 90% width | Mean measured wall time (s) |
|---|---:|---:|---:|---:|---:|---:|
| fixed_direct | 3 | 1.238 | 8.626 | — | — | 732.0 |
| active_solvai_bq | 3 | 1.888 | 9.061 | 0.875 | 7.722 | 732.0 |
| fixed_solvai_bq | 3 | 1.888 | 9.061 | 0.875 | 7.722 | 732.0 |
| oracle_non_deployable | 3 | 1.888 | 9.061 | 0.875 | 7.722 | 732.0 |
| random_solvai_bq | 3 | 1.888 | 9.061 | 0.875 | 7.722 | 732.0 |
| uniform_solvai_bq | 3 | 1.888 | 9.061 | 0.875 | 7.722 | 732.0 |
| generic_bq | 3 | 1.975 | 7.310 | 0.875 | 8.037 | 732.0 |
| oracle_non_deployable | 5 | 0.337 | 8.103 | 1.000 | 6.659 | 835.5 |
| active_solvai_bq | 5 | 1.701 | 10.093 | 0.750 | 5.406 | 829.0 |
| uniform_solvai_bq | 5 | 1.701 | 10.093 | 0.750 | 5.406 | 829.0 |
| generic_bq | 5 | 1.710 | 7.506 | 0.750 | 4.421 | 829.0 |
| fixed_direct | 5 | 1.727 | 4.947 | — | — | 828.7 |
| random_solvai_bq | 5 | 1.845 | 8.826 | 0.850 | 6.831 | 829.9 |
| fixed_solvai_bq | 5 | 3.044 | 6.226 | 0.750 | 7.586 | 828.7 |
| oracle_non_deployable | 7 | 0.068 | 9.274 | 1.000 | 6.007 | 931.6 |
| generic_bq | 7 | 1.114 | 7.896 | 0.750 | 3.732 | 927.3 |
| active_solvai_bq | 7 | 1.608 | 11.259 | 0.875 | 4.893 | 926.7 |
| uniform_solvai_bq | 7 | 1.625 | 9.510 | 0.750 | 5.167 | 927.5 |
| fixed_direct | 7 | 1.691 | 5.903 | — | — | 925.6 |
| random_solvai_bq | 7 | 1.956 | 8.411 | 0.744 | 6.024 | 927.3 |
| fixed_solvai_bq | 7 | 2.548 | 6.873 | 0.500 | 5.202 | 925.6 |
| oracle_non_deployable | 9 | 0.088 | 9.979 | 1.000 | 5.130 | 1026.7 |
| generic_bq | 9 | 0.911 | 9.523 | 0.750 | 2.971 | 1026.6 |
| uniform_solvai_bq | 9 | 1.279 | 11.614 | 0.875 | 4.571 | 1027.3 |
| active_solvai_bq | 9 | 1.438 | 13.345 | 0.750 | 4.331 | 1026.5 |
| random_solvai_bq | 9 | 1.738 | 8.086 | 0.744 | 5.212 | 1025.2 |
| fixed_solvai_bq | 9 | 1.799 | 5.942 | 0.625 | 5.074 | 1020.2 |
| fixed_direct | 9 | 1.902 | 6.495 | — | — | 1020.2 |
| fixed_direct | 15 | 0.000 | 0.000 | — | — | 1316.0 |
| generic_bq | 15 | 0.475 | 0.000 | 0.875 | 1.925 | 1316.0 |
| active_solvai_bq | 15 | 1.214 | 0.000 | 0.750 | 3.382 | 1316.0 |
| fixed_solvai_bq | 15 | 1.214 | 0.000 | 0.750 | 3.382 | 1316.0 |
| oracle_non_deployable | 15 | 1.214 | 0.000 | 0.750 | 3.382 | 1316.0 |
| random_solvai_bq | 15 | 1.214 | 0.000 | 0.750 | 3.382 | 1316.0 |
| uniform_solvai_bq | 15 | 1.214 | 0.000 | 0.750 | 3.382 | 1316.0 |

## Interpretation

The frozen positive criterion was not met. Direction B is killed for this protocol and Direction C is not launched.

The active policy first reached 0.20 kcal mol⁻¹ MAE at no tested windows. Comparator minima were `{"curvature_solvai_bq": null, "fixed_direct": 15, "fixed_solvai_bq": null, "generic_bq": null, "random_solvai_bq": null, "uniform_direct": 15, "uniform_solvai_bq": null}`.

The oracle is non-deployable and cannot support a method claim. All twelve molecules were already development-exposed in the parent project; only their twelve added dense-window responses were prospective.
