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

