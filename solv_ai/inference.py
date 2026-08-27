"""Public SMILES-only inference API for the frozen SolvAI model stack."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np

from .features import canonicalize, descriptor_frame
from .teachers import response_features

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL_DIR = ROOT / "models/final"


def predict_smiles(
    smiles: list[str], model_dir: str | Path = DEFAULT_MODEL_DIR
) -> tuple[np.ndarray, np.ndarray]:
    """Return hydration free energy and seed-ensemble spread in kcal/mol.

    The only molecular input is SMILES. All response priors are predictions from
    bundled structure models; no molecular simulation or experimental lookup is
    performed.
    """
    if not smiles:
        raise ValueError("At least one SMILES string is required")
    model_dir = Path(model_dir)
    canonical = [canonicalize(value)[0] for value in smiles]
    descriptors = descriptor_frame(canonical)
    physical = response_features(canonical, descriptors, model_dir)
    artifact: dict[str, Any] = joblib.load(model_dir / "head.joblib")
    inputs = np.column_stack(
        [descriptors[artifact["descriptor_columns"]].to_numpy(dtype=float), physical]
    )
    expected = len(artifact["descriptor_columns"]) + len(artifact["teacher_columns"])
    if inputs.shape[1] != expected:
        raise AssertionError(f"Inference schema has {inputs.shape[1]} columns; expected {expected}")
    members = np.column_stack([model.predict(inputs) for model in artifact["models"]])
    return members.mean(axis=1), members.std(axis=1)
