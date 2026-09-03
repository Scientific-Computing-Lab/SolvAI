# Phase 0 preservation, reproduction and response inventory

## Parent state

- Repository: `/home/galoren/Freecurve_AI_Solvation/release/SolvAI`
- Starting branch: `main`
- Starting commit: `531f6cfd21e319c951b461c9ef24fa754790f91d`
- Starting working tree: clean
- Active branch: `active-solvai`
- Parent commit discrepancy: none

The parent paper, Supplement, freezes, model artifacts and results remain unmodified.

## Canonical parent locations

| Role | Path |
|---|---|
| Main manuscript | `paper/main.tex`, `paper/main.pdf` |
| Supplement | `paper/supplementary/supplementary.tex`, `paper/supplementary/supplementary.pdf` |
| Confirmatory report | `reports/CONFIRMATORY_ANALYSIS.md` |
| Tier-A report | `reports/TIER_A_EXTERNAL_VALIDATION.md` |
| Provenance | `repro/DATA_PROVENANCE.md` |
| Parent model | `models/final/` |
| Parent fixed predictions | `results/confirmatory/standardized_exclusion_endpoint_predictions.parquet` |
| Parent molecule predictions | `results/oof_predictions.parquet` |
| ARROW benchmark | `data/benchmark/arrow_solvation_master.parquet` |
| Historical PIMD2 response tables | workspace `results/pimd2_multilambda_teacher*.parquet` |
| Raw PIMD2 probe outputs | workspace `results/physics_probes/pimd2_lambda*_5ps/` |

## Frozen parent hashes

| Artifact | SHA-256 |
|---|---|
| Main PDF | `7bf883c62828e08f25b65fded21a4a149dc776f198f0a26ee03394b6fd52f041` |
| Supplement PDF | `58b46b69f5b6d23a21e0eece00103b9c904f4227874aad345a6bcdc87cc81d89` |
| Parent release manifest | `645da9d30ece89d59b2bde47ffbd0a216b53754ea00ea69f93e2b42975c8d1ad` |
| Parent model manifest | `d64d6cbca36b9b5cc672c9b7622c5b72dc1e9409cfd0daaf6a7a6d36c5ec64c7` |
| Parent paper metrics | `cb8c788234ded4764aad85900fc1e82c7c1840d8140f5c6c7847419809e3d946` |

All nine files named by the model manifest reproduced their recorded size and SHA-256.

## Parent reproduction

The two packaged SMILES-only smoke predictions reproduced within 2e-6 kcal mol-1. The corrected fixed matched values were independently recomputed from the canonical prediction table:

- matched structure-only: 0.3033484656 kcal mol-1;
- frozen SolvAI: 0.2022340672 kcal mol-1;
- N = 85.

Machine-readable output: `results/phase0/parent_reproduction.json`.

## Host and software

- CPU: Intel Core i7-4930K, 6 physical / 12 logical cores.
- RAM: 62 GiB.
- GPU: NVIDIA GeForce RTX 3090, 24,576 MiB, driver 535.309.01.
- Disk at start: 1.8 TiB filesystem, 116 GiB free (94% used).
- Parent Python: 3.11.15; NumPy 2.4.6; pandas 3.0.5; pyarrow 25.0.1; SciPy 1.17.1; scikit-learn 1.7.2; RDKit 2026.03.5; PyTorch 2.5.1+cu124; Chemprop 2.2.4; LightGBM 4.7.0.
- `gpytorch` and OpenMM are not installed in the parent environment. The first gate does not require either.

Exact tool paths and complete system output are stored in `data/manifests/environment.json`.

## Existing actual response observations

Three private local PIMD2 campaigns use 5 ps production at nominal lambda 0.1, 0.5 and 0.9:

| Campaign | Attempted | Successful |
|---|---:|---:|
| lambda 0.1 | 76 | 72 |
| lambda 0.5 | 76 | 74 |
| lambda 0.9 | 76 | 73 |

Exactly 72 ARROW molecules have complete three-point observations. All 76 attempted probe identities map to the 85-row water benchmark; no identity is unresolved. The complete cohort spans 14 diagnostic families.

The inventory hashes 1,356 configurations, logs, timings and energy files (25,729,358 bytes). It parses 35,680 sequential-prefix response summaries at 10%, 20%, 40%, 70% and 100% of production, with naive and five-contiguous-block uncertainty estimates. Raw paths remain outside the release tree.

## Dense-response availability

No compatible dense multi-window **population** was located. One 11-window toluene software-test dataset and one legacy 21-window directory were found, but neither can support population replay: the first has one molecule, and the second lacks resolved molecule/protocol provenance. Consequently:

- Phase 1 endpoint and trajectory-prefix gates are feasible with zero new simulation.
- Same-Hamiltonian hidden-window population reconstruction and Phase 2 adaptive replay are not currently feasible from existing local data alone.
- If Phase 1 passes, the next step is a predeclared small dense sentinel pool, not reinterpretation of unrelated public alchemical data.

## Files produced

- `data/identity/probe_identity_manifest.{csv,parquet}`
- `data/manifests/response_case_inventory.csv`
- `data/manifests/response_file_manifest.csv`
- `data/manifests/probe_split_assignments.csv`
- `data/manifests/dense_response_candidates.json`
- `results/phase0/response_prefix_blocks.parquet`
- `results/phase0/response_inventory_summary.json`

These outputs distinguish raw files, historical derived tables and new deterministic summaries.

