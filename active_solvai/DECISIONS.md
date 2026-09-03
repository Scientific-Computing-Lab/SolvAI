# Decision log

This file is append-only. Amendments cite the affected prospective freeze and state what results were already known.

## 2026-09-03 — Isolate Active SolvAI

- **Decision:** create branch `active-solvai` and keep all new artifacts below `active_solvai/`.
- **Reason:** preserve the completed SolvAI release at parent commit `531f6cfd21e319c951b461c9ef24fa754790f91d`.
- **Alternatives rejected:** modifying the parent manuscript/results tree or merging into `main` before approval.
- **Evidence available:** parent repository clean at the expected commit; blueprint hash recorded in the ingestion report.

## 2026-09-03 — Start with existing PIMD2 observations

- **Decision:** use the existing three 5 ps PIMD2 campaigns at nominal lambda 0.1, 0.5 and 0.9 for the first gate; launch no new simulation before that gate is frozen and scored.
- **Reason:** the blueprint requires the cheapest decisive test first, and raw logs plus `.ene` time series are present locally.
- **Known limitation:** these are sparse, short, protocol-specific observations, not dense PIMD8 response curves.

## 2026-09-03 — Kill endpoint correction; permit one dense sentinel

- **Evidence visible:** the frozen Phase 1 endpoint result had been scored. The three-point actual-minus-predicted response increased five-repeat MAE by 0.003482 kcal mol⁻¹ (95% paired interval +0.000960 to +0.006093) and did not beat shuffled residuals.
- **Decision A:** kill empirical endpoint residual correction for this 5 ps PIMD2 protocol. No λ subset, component block or chemical family may rescue it.
- **Separate reconstruction evidence:** molecule-conditioned Gaussian interpolation reduced held-out response-point error versus a generic Gaussian for several predeclared subsets with approximately nominal coverage, although absolute errors and intervals remained large. No compatible dense response population exists locally.
- **Decision B:** permit one bounded, prospectively frozen dense PIMD2 sentinel acquisition to test same-Hamiltonian reconstruction. This does not reopen the endpoint gate.
- **Decision C:** do not launch adaptive multi-fidelity escalation unless the dense sentinel first establishes useful reconstruction against fixed and generic comparators.
