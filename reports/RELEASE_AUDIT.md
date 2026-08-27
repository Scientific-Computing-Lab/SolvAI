# SolvAI publication-release audit

## Scientific freeze

All manuscript values are recomputed by `solv_ai.paper_metrics` from canonical
molecule-level outputs. Hard assertions protect the reference-set size, core
MAEs, external-source counts and model feature schema. The authoritative files
are `results/paper_metrics.json`, `results/paper_metrics.csv` and
`reports/PAPER_FREEZE.md`.

The release distinguishes three quantities that must not be conflated:

- fixed strict five-fold OOF: 0.19705 kcal/mol;
- nested-selection five-fold OOF: 0.19931 kcal/mol;
- five independent fixed-model repeats: 0.20375 ± 0.00493 kcal/mol.

The claim is PIMD8-level accuracy on the ARROW 85-solute reference set, not
robust universal sub-0.20 performance.

## Leakage and inference boundary

The final independent audit passed. Every supervised external source is
disjoint from all 85 reference connectivities. Redistributed tables are
recanonicalized at release verification; non-redistributed CombiSolv-QM and
mixed public-hydration tables are represented by their frozen clean-room audit,
exact SHA-256 values and source URLs. The inference schema contains 2,265
structure descriptors and 15 predicted response priors. No experimental target,
ARROW or PIMD value, trajectory, probe, family, scaffold or fold field is
present.

## Artifact and runtime

All bundled model hashes and two SMILES-only predictions reproduce within
2 × 10⁻⁶ kcal/mol. Median warm latency on the release workstation is 17.522 s
for one molecule; batch-32 throughput is 1.79 molecules/s. The RTX 3090 was
present but not used by packaged CPU inference. Published PIMD8 work estimates
are reported separately because they were measured on different hardware.

## Publication and reproducibility

- Main manuscript and Supplementary Information compile from source.
- Five main and eight Extended Data figures are regenerated from frozen files.
- The master experiment table preserves non-improving routes.
- Four scientific review perspectives were completed, followed by one
  post-revision pass.
- Automated claim and credential scans pass.
- Human-only declarations remain isolated in `submission/AUTHOR_TODOS.md`.

## Redistribution

Large original archives are not copied. Processed CC BY/CC BY-SA tables are
included with attribution and Git LFS; CombiSolv-QM and merged endpoint-label
rows are excluded where no standalone data licence is stated. The trained
checkpoints, source DOIs, exact processed hashes and opt-in download commands
remain available.
