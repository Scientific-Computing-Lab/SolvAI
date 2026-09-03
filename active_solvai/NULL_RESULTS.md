# Null-result register

This file is append-only. Each entry identifies the hypothesis, frozen configuration, run IDs, exact result, and scientific interpretation.

## Historical context (not an Active SolvAI result)

The completed SolvAI campaign found that **structure-predicted** PIMD2 response features did not improve its matched hydration endpoint. That result does not test whether an **actual observation from the target molecule**, or its residual relative to a structure-predicted prior, contains useful information.



## AS-P1-GATE-001 — actual PIMD2 observation gate (2026-09-03)

- **Frozen before scoring:** commit `a0dd986`.
- **Question:** do three actual 5 ps PIMD2 SYSTEM dH/dλ observations, especially actual-minus-structure-predicted residuals, improve frozen SolvAI?
- **Result:** no. Five-repeat mean MAE changed from 0.186374 to 0.189856 kcal mol⁻¹; paired difference +0.003482 (95% CI +0.000960, +0.006093). The aligned residual was not better than shuffled residuals.
- **Mechanistic diagnostic:** structure-conditioned interpolation reduced some hidden-point errors, but errors remained ≥2.4 kcal mol⁻¹ with broad intervals; three points do not constitute a dense reconstruction benchmark.
- **Decision:** kill Direction A for this probe protocol. Hold Direction C. The limited hidden-point signal permits one separately frozen dense sentinel test of Direction B; it cannot alter this endpoint null.

## AS-P2-DENSE-PIMD2-001 — prospective dense sentinel (2026-09-03)

- **Frozen before prospective simulation:** scientific protocol commit
  `512aad3`; leakage amendment `c981630`; calibration lock `3e69d10`.
- **Question:** can a structure-conditioned curve prior and Bayesian acquisition
  reconstruct a short-window, same-Hamiltonian dense PIMD2 integral more
  accurately than cost-matched fixed, uniform, random, curvature-only and
  generic-BQ policies?
- **Prospective evidence:** 8 sentinel molecules, 15 λ values each, with 12 new
  5-ps PIMD2 windows per molecule. All 96 prospective windows passed QC on the
  first attempt.
- **Result:** no. At five windows, active BQ MAE was 1.701 kcal mol⁻¹ versus
  1.153 for the strongest comparator (uniform direct), a paired difference of
  +0.548 (90% CI −0.101, +1.159), with 2/8 molecules improved. At seven windows,
  active BQ was 1.608 versus 1.092, a paired difference of +0.516 (90% CI
  −0.392, +1.491), with 3/8 improved.
- **Calibration:** 90% coverage was 0.750 and 0.875 at five and seven windows,
  but interval widths remained 5.406 and 4.893 kcal mol⁻¹. No Bayesian method
  reached the frozen 0.10-kcal mol⁻¹ half-width stopping threshold.
- **Diagnostic:** the non-deployable same-curve oracle reached 0.337 and 0.068
  kcal mol⁻¹ at five and seven windows. A later held-out-block audit reversed
  this apparent advantage and found unstable selected sets; stable
  molecule-specific placement headroom is therefore not established.
- **Decision:** kill Direction B for this protocol. Direction C is not launched
  because it was prospectively contingent on Direction B passing. Tier-B
  remains unopened.

## AS-V2-REPLICA-POWER-001 — prospective simulation power gate (2026-09-03)

- **Frozen before calculation:** commit `52827f851ed0d1a540b0492c63c3df269cf774e1`.
- **Hypothesis:** two independent 50-ps replicas across eight molecules and all 15 windows can resolve a stable mean oracle-placement advantage of 0.30 kcal mol$^{-1}$ at five and seven windows.
- **Result:** no. Power was 0.00000 and 0.00001 at five and seven windows; the dense-reference reliability probability was 0.44602. The registered launch thresholds were 0.80 for all three quantities.
- **Compute:** zero new simulation, zero bead-windows and zero GPU-hours. The proposed 36-GPU-hour campaign was stopped prospectively.
- **Interpretation:** observed molecule-level heterogeneity, not only short-trajectory noise, prevents the small design from establishing reproducible placement headroom. No experimental A/B/C decision is assigned because the gated campaign did not run.
