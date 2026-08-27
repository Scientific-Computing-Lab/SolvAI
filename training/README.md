# Frozen training pipeline

This directory preserves the exact Python sources used for the final teacher
and endpoint stages. It is intentionally outside the import path of the compact
inference package.

Install the full environment with:

```bash
uv sync --extra dev --extra training
export PYTHONPATH="$PWD/training"
```

The scripts expect the directory layout encoded in `arrow_distill/data.py`.
Downloaded Zenodo files live under `data/external/zenodo/<record>/`; archives
must be extracted in place. Publisher supplements are placed at the explicit
paths named in `scripts/build_global_data_catalog.py`. The immutable source and
processed hashes are in `../data/manifests/training_source_manifest.json`.
The mixed experimental endpoint-label tables are reconstructed from the frozen
FreeSolv, CombiSolv-Exp and SoluteML sources recorded in
`../data/manifests/endpoint_label_manifest.json`; they are deliberately not
redistributed because one component has no standalone data licence.

The core final sequence was:

```bash
python training/scripts/build_global_data_catalog.py
python training/scripts/prepare_soluteml_hydration.py
python training/scripts/prepare_openff_alchemical_teacher.py
python training/scripts/prepare_implicit_solvent_teacher.py
python training/scripts/prepare_molsolv_smd_teacher.py
python training/scripts/pretrain_combisolv_qm.py
python training/scripts/pretrain_molsolv_smd.py
python training/scripts/prepare_confsolv_water_teacher.py
python training/scripts/train_confsolv_water_teacher.py
python training/scripts/serialize_tree_teacher_artifacts.py
python training/scripts/run_nested_smd_teacher_confirmation.py --no-static-arrow --confsolv-response
python training/scripts/train_final_structure_only_artifact.py
```

Teacher pretraining is expensive and source-specific; see
`../repro/FULL_REPRODUCTION.md` for expected resources. Frozen held-out
predictions, rather than a rerun of this pipeline, are the canonical basis for
the publication metrics.

Private Freecurve simulator source and trajectories are not included. They are
not needed to rebuild or run the selected SolvAI artifact; no PIMD8 feature was
retained in that model.
