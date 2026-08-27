"""Refit deployment teacher models using already selected CV hyperparameters.

The expensive cross-validation and leakage-safe prediction tables already exist.
This script only reproduces each final full-source structure surrogate for the
standalone SMILES predictor; it does not perform model selection or evaluation.
"""

from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from sklearn.ensemble import ExtraTreesRegressor


def estimator(seed: int, leaf: int, trees: int) -> ExtraTreesRegressor:
    return ExtraTreesRegressor(
        n_estimators=trees,
        max_features=0.7,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )


def save_artifact(
    directory: str,
    source: pd.DataFrame,
    source_x: np.ndarray,
    descriptor_columns: list[str],
    targets: list[str],
    metadata_targets: dict[str, dict[str, object]],
    trees: int,
) -> None:
    models = {}
    for target_index, target in enumerate(targets):
        labeled = source[target].notna().to_numpy()
        leaf = int(metadata_targets[target]["selected_min_samples_leaf"])
        model = estimator(20260826 + 1009 * target_index, leaf, trees)
        model.fit(source_x[labeled], source.loc[labeled, target].to_numpy(dtype=float))
        models[target] = model
        print(directory, target, int(labeled.sum()), leaf, flush=True)
    output = ROOT / "models" / directory
    output.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "models": models,
            "descriptor_columns": descriptor_columns,
            "target_columns": targets,
        },
        output / "models.joblib",
        compress=3,
    )


def main() -> None:
    processed = ROOT / "data/processed"
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    public_features = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    descriptor_columns = [
        column for column in public_features if column.startswith(("rdkit__", "morgan2__"))
    ]
    x_public = public_features[descriptor_columns].to_numpy(dtype=np.float32)
    public_position = {
        key: index for index, key in enumerate(public.inchi_connectivity_key.astype(str))
    }

    abraham = pd.read_parquet(processed / "soluteml_abraham_nonbenchmark.parquet")
    abraham_features = pd.read_parquet(processed / "soluteml_abraham_rdkit_morgan_features.parquet")
    if not abraham.molecule_id.equals(abraham_features.molecule_id):
        raise AssertionError("Abraham source feature order changed")
    abraham_targets = [f"abraham_{name}" for name in "esabl"]
    abraham_metadata = json.loads((processed / "soluteml_hydration_metadata.json").read_text())[
        "abraham_targets"
    ]
    save_artifact(
        "soluteml_abraham_teacher",
        abraham,
        abraham_features[descriptor_columns].to_numpy(dtype=np.float32),
        descriptor_columns,
        abraham_targets,
        abraham_metadata,
        160,
    )

    for directory, source_file, metadata_file, targets, trees in (
        (
            "openff_alchemical_teacher",
            "openff_alchemical_nonbenchmark.parquet",
            "openff_alchemical_teacher_metadata.json",
            [
                "openff23_dg",
                "openff23_solvent_decoupling",
                "openff23_vacuum_decoupling",
                "openff23_exp_residual",
                "openff23_minus_legacy_gaff",
                "ash_consensus_dg",
                "ash_method_spread",
            ],
            500,
        ),
        (
            "implicit_solvent_teacher",
            "implicit_solvent_nonbenchmark.parquet",
            "implicit_solvent_teacher_metadata.json",
            [
                "gbn2_alchemical_dg",
                "neural_implicit_dg",
                "gbn2_exp_residual",
                "neural_implicit_exp_residual",
                "neural_minus_gbn2",
            ],
            500,
        ),
    ):
        source = pd.read_parquet(processed / source_file)
        positions = np.asarray(
            [public_position[str(key)] for key in source.connectivity_key.astype(str)],
            dtype=int,
        )
        metadata = json.loads((processed / metadata_file).read_text())["targets"]
        save_artifact(
            directory,
            source,
            x_public[positions],
            descriptor_columns,
            targets,
            metadata,
            trees,
        )


if __name__ == "__main__":
    main()
