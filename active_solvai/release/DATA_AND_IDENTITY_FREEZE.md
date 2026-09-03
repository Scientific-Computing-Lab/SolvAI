# Active SolvAI data and identity freeze

**Date/timezone:** 2026-09-03 UTC  
**Parent SolvAI commit:** `531f6cfd21e319c951b461c9ef24fa754790f91d`  
**Active branch:** `active-solvai`  
**Result visibility:** no Active SolvAI performance result had been scored when this document was written.

## Scientific scope

The first gate uses neutral water-solvation molecules from ARROW-85 that have successful actual PIMD2 observations at every nominal lambda in {0.1, 0.5, 0.9}. Eligibility is determined only by parent identity and simulation completion, never endpoint error.

## Identity

The stable molecule key is the parent full InChIKey (`molecule_id`). Canonical isomeric SMILES, connectivity block, diagnostic family and scaffold are inherited from `data/benchmark/arrow_solvation_master.parquet`. All attempted probe names map uniquely to those identities. No molecule is removed because of prediction error, family or response magnitude.

## Frozen cohort

- ARROW water reference: 85 molecules.
- Probe attempts: 76 unique molecules.
- Successful lambda 0.1: 72.
- Successful lambda 0.5: 74.
- Successful lambda 0.9: 73.
- Complete matched primary cohort: 72.

The primary panel uses these same 72 molecules for every method. Partial observations from the remaining four attempted molecules may appear only in explicitly labelled missingness diagnostics; they cannot enlarge a favorable method row.

## Splits

- Primary fixed folds: parent `fold_random` assignment.
- Repeat folds: exact parent assignments in `data/manifests/probe_split_assignments.csv` for seeds 314159, 271828, 161803, 141421 and 173205.
- Chemical checks: inherited `fold_family` and `fold_scaffold`, subject to the Phase 1 freeze’s minimum-training-size rule.

No split may be altered after result visibility.

## Raw and derived data

| File | SHA-256 |
|---|---|
| workspace `results/pimd2_multilambda_teacher.parquet` | `56b9607a6e53d22c7ec1d504102594141314533e4ca9b90277cc5516675998de` |
| workspace `results/pimd2_multilambda_teacher_long.parquet` | `6e6d445ab8e0ece61aeaa519d73b537f5a077ddda70dd87660f785a2e14e0fef` |
| `data/identity/probe_identity_manifest.parquet` | `4fd272eaa4b59e18ef31a9ae094fd37bfb5d0232898dd0846c4a4fca318e3c11` |
| `data/manifests/response_case_inventory.csv` | `f712ce936835d52e12b7863ceaa42458509ebedd5b629dc1dab9de903fc23c3f` |
| `data/manifests/response_file_manifest.csv` | `f7d016937816d4b5c27b8591dc5305c17c7eaf404d6799ea3bf754ff6b029fdc` |
| `data/manifests/probe_split_assignments.csv` | `f2fb6e9334f7f95da5fd5a52c0d1c1a377cad12fe49528157bc96ddc306fe0c8` |
| `results/phase0/response_prefix_blocks.parquet` | `4d316e5695b49930451e6eb5c2073d26cc87f163e5eac04552dbea4517e8fbae` |

The file manifest stores absolute authorized-local paths, sizes and hashes for 1,356 source files. Raw simulation outputs are not copied into Git.

## Dense replay constraint

The local inventory does not contain a compatible dense multi-window population. Phase 1 therefore tests experimental endpoint value and trajectory-prefix stability. Dense reconstruction is not silently approximated from the three points. A compatible dense sentinel dataset must be prospectively frozen or obtained from collaborators before Direction B/Phase 2 can claim population reconstruction.

## Amendment rule

This file is immutable after commit. Any change is a separate dated amendment that identifies the information already observed, exact changed field and scientific consequence.

