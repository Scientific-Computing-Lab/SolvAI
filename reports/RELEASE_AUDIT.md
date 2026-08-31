# SolvAI publication-release audit

## Scientific freeze

The preregistered confirmatory package supersedes the exploratory campaign headline.
All manuscript values are recomputed by `solv_ai.paper_metrics` from molecule-level
outputs. The authoritative files are `results/paper_metrics.json`,
`results/paper_metrics.csv`, `reports/PAPER_FREEZE.md` and
`reports/CONFIRMATORY_ANALYSIS.md`.

The matched primary comparison is 0.30335 kcal/mol for the descriptor-only endpoint
and 0.20223 kcal/mol for SolvAI. The paired change is −0.10111 kcal/mol (95% interval,
−0.2151 to −0.0203). Five complete partitions give 0.31295 ± 0.00421 and
0.20737 ± 0.00444 kcal/mol, respectively. The reconstructed ARROW/PIMD8 comparator is
0.20484 kcal/mol. The claim is PIMD8-level accuracy on this reference chemistry, not
superiority to PIMD8 or broadly sub-0.20 performance.

A prospectively frozen Sander Tier-A analysis retains 220 endpoint-disjoint molecules,
including 97 also absent from all response-teacher sources. The matched comparison is
1.53165 versus 1.15255 kcal/mol on N=220 and 2.13830 versus 1.53560 kcal/mol on strict
N=97; paired 95% intervals exclude zero in both. This supports transfer of the
response-prior advantage while showing substantially weaker absolute accuracy on
broader chemistry.

## Leakage and inference boundary

The confirmatory identity audit adds fragment-parent, uncharging and tautomer
standardization to exact connectivity. It found and removed 2 CombiSolv-QM, 32
MolSolv and 22 ConfSolv benchmark equivalents, then refitted the affected teachers
without changing the split membership of any retained source molecule. The 1,280
external endpoint labels contain no exact or standardized benchmark identity.

The deployed schema contains 2,265 structure features and 15 structure-predicted
response priors. No experimental target, ARROW or PIMD value, trajectory, probe,
family, scaffold or fold field is present. PIMD-derived candidate supervision was not
retained.

## Artifact and runtime

All bundled hashes and two SMILES-only predictions reproduce within 2 × 10⁻⁶
kcal/mol. Median warm latency is 15.294 s for one molecule; batch-32 throughput is
2.022 molecules/s (0.494 s per molecule) on an Intel i7-4930K CPU. The installed RTX
3090 was not used. Published PIMD8 total-work estimates are reported separately
because they were measured on different hardware.

## Publication and reproducibility

- Four main and six Extended Data figures are regenerated from frozen files.
- Supplementary Information contains methods, notes and tables but no figures.
- Machine-readable Supplementary Data preserve predictions, splits and all campaign
  experiments.
- Human-only declarations remain isolated in `submission/AUTHOR_TODOS.md`.

## Redistribution

Large source archives are not copied. Processed redistributable tables are included
with attribution and Git LFS; restricted publisher supplements are represented by
source links, hashes, filtering logic and opt-in download commands.
