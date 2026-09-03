# Null-result register

This file is append-only. Each entry must identify the hypothesis, frozen configuration, run IDs, exact result, and scientific interpretation. No Active SolvAI result has been scored yet.

## Historical context (not an Active SolvAI result)

The completed SolvAI campaign found that **structure-predicted** PIMD2 response features did not improve its matched hydration endpoint. That result does not test whether an **actual observation from the target molecule**, or its residual relative to a structure-predicted prior, contains useful information.



## AS-P1-GATE-001 — actual PIMD2 observation gate (2026-09-03)

- **Frozen before scoring:** commit `a0dd986`.
- **Question:** do three actual 5 ps PIMD2 SYSTEM dH/dλ observations, especially actual-minus-structure-predicted residuals, improve frozen SolvAI?
- **Result:** no. Five-repeat mean MAE changed from 0.186374 to 0.189856 kcal mol⁻¹; paired difference +0.003482 (95% CI +0.000960, +0.006093). The aligned residual was not better than shuffled residuals.
- **Mechanistic diagnostic:** structure-conditioned interpolation reduced some hidden-point errors, but errors remained ≥2.4 kcal mol⁻¹ with broad intervals; three points do not constitute a dense reconstruction benchmark.
- **Decision:** kill Direction A for this probe protocol. Hold Direction C. The limited hidden-point signal permits one separately frozen dense sentinel test of Direction B; it cannot alter this endpoint null.
