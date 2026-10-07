# Source-level residual construction

These two unmodified archived scripts document the original OpenFF/GBn2 target
definitions and public-row feature construction. They form each experimental
residual against the original physical calculation. Public molecules overlapping
the relevant source receive source-level five-fold OOF predictions; other public
and benchmark rows receive full-source predictions. The same source folds select
minimum leaf size, so their OOF errors are selection scores, not nested unbiased
validation estimates. These folds are not rebuilt around later endpoint splits.

Original SHA-256:

- `prepare_openff_alchemical_teacher.py`: `75c41781468a9f980b5cc0de2b983137052e18a37e3610c409723eaf39f19a73`
- `prepare_implicit_solvent_teacher.py`: `c4125ea1444f2b629ff65c9b7e4755e92864c799337ef338a2f00f0c2cb00318`

They are implementation records, not standalone training entry points. They use
the historical `arrow_distill.data` adapter for workspace root, the FreeSolv JSON
path and canonicalization, plus acquired original source files and precomputed
feature tables. Adapt these inputs explicitly when reconstructing source training;
do not point them at an existing release directory. Acquisition and redistribution
terms are in `../DATA_PROVENANCE.md`; source metadata are under `data/manifests/`.
No original source labels are added by this archive. Current public inference
continues to use `solv_ai/teachers.py` and the frozen release weights.
