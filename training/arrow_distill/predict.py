"""Inference interface for the zero-simulation ARROW distillation ensemble."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from .data import ROOT, canonicalize, rdkit_descriptor_frame
from .sklearn_compat import repair_legacy_simple_imputers

DEFAULT_ARTIFACT = ROOT / "models/final_ensemble/final_ensemble.joblib"


def load_artifact(path: str | Path = DEFAULT_ARTIFACT) -> dict[str, Any]:
    return repair_legacy_simple_imputers(joblib.load(path))


def predict_smiles(
    smiles: list[str],
    artifact: dict[str, Any] | None = None,
    static_features: list[dict[str, float]] | None = None,
) -> np.ndarray:
    """Predict hydration free energies in kcal/mol without simulation.

    Static ARROW aggregates are optional. When absent, the model follows its
    learned missing-static path. No trajectory, ARROW energy, or PIMD result is
    consumed.
    """
    if artifact is None:
        artifact = load_artifact()
    canonical = [canonicalize(item)[0] for item in smiles]
    records = pd.DataFrame(
        {
            "molecule_id": [f"query_{index}" for index in range(len(smiles))],
            "canonical_smiles": canonical,
        }
    )
    features = rdkit_descriptor_frame(records)
    x = features[artifact["feature_columns"]].to_numpy(dtype=float)
    static = np.full((len(smiles), len(artifact["static_columns"])), np.nan)
    static[:, 0] = 0.0
    if static_features is not None:
        if len(static_features) != len(smiles):
            raise ValueError("static_features and smiles must have equal lengths")
        for row_index, values in enumerate(static_features):
            for column_index, column in enumerate(artifact["static_columns"]):
                if column in values:
                    static[row_index, column_index] = float(values[column])
    x_static = np.hstack([x, static])
    component_predictions = np.column_stack(
        [
            artifact["components"]["public"].predict(x),
            artifact["components"]["static"].predict(x_static),
            artifact["components"]["distilled"].predict(x),
        ]
    )
    return component_predictions @ np.asarray(artifact["weights"])
