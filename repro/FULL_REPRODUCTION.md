# Full reproduction

The public release separates fast artifact reproduction from training-time
physics reconstruction. The final inference stack is fully included; rebuilding
every teacher requires large external data, substantial storage and a GPU.

## Scope

The full workflow comprises:

1. download each source at the version in `DATA_PROVENANCE.md`;
2. validate its published checksum and licence;
3. canonicalize solutes with RDKit;
4. remove exact and standardized ARROW-85 equivalents before supervised fitting;
5. build the MolSolv, ConfSolv, CombiSolv-QM, Abraham, OpenFF and implicit
   response teachers;
6. preserve each retained molecule's frozen teacher split membership and refit the
   affected response teachers;
7. generate response predictions and fit the preregistered confirmatory endpoints;
8. refit the frozen deployment artifact only after evaluation.

The exact final-stage training sources are preserved under `training/`; the full
40-repository Freecurve snapshot and unrelated exploratory scripts are not
redistributed because they include private code and third-party material. The
release preserves the executable final inference implementation, training
configuration, filtered redistributable tables, source hashes, molecule-level
held-out predictions and every manuscript result.

After arranging external archives as described in `training/README.md`, use:

```bash
uv sync --extra dev --extra training
export PYTHONPATH="$PWD/training"
```

## Expected resources

- CPU: at least 16 cores recommended for RDKit preprocessing and tree teachers.
- GPU: one CUDA-capable accelerator with 20 GB or more memory for D-MPNN
  teachers.
- Storage: at least 200 GB for source archives and intermediates.
- Time: hours to days depending on download rate and accelerator.

No new PIMD calculation is needed to rebuild the selected model. The PIMD2 and
PIMD8 results in the supplement are candidate-supervision analyses and are not
part of the final artifact.

## Non-redistributed inputs

CombiSolv supplementary data have no standalone data licence in the archive and
are therefore not copied into this repository. This applies to the
CombiSolv-QM response table and to merged public hydration tables containing
CombiSolv-Exp records. `scripts/download_data.py` retrieves these sources only
after the researcher accepts the source terms. The manifests record the exact
filtered-file hashes used in this study. To reconstruct the merged endpoint
labels, download SoluteML, FreeSolv and CombiSolv-Exp, then run
`training/scripts/prepare_soluteml_hydration.py`. The final trained checkpoints
are included for quick reproduction.

## Frozen evaluation

The definitive inputs are the molecule-level files under `results/confirmatory/`,
not predictions from the all-data deployment refit. Historical campaign files under
`results/predictions`, `results/robustness` and `results/ablations` are retained for
provenance and are labelled as such. `solv_ai.paper_metrics` is the only module
authorized to write manuscript metrics.
