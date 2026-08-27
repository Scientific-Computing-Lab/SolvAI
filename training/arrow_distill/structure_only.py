"""SMILES-only inference for the frozen physics-distilled hydration model."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from .data import ROOT, canonicalize, rdkit_descriptor_frame

DEFAULT_MODEL_DIR = ROOT / "models/final_structure_only"


def _chemprop_prediction(
    smiles: list[str],
    checkpoint: Path,
    target_prefix: str,
    *,
    rdkit_2d: bool,
) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="arrow_structure_only_") as directory:
        root = Path(directory)
        input_path = root / "input.csv"
        output_path = root / "output.csv"
        pd.DataFrame({"smiles": smiles}).to_csv(input_path, index=False)
        executable = ROOT / ".venv/bin/chemprop"
        command = [
            str(executable),
            "predict",
            "-i",
            str(input_path),
            "-s",
            "smiles",
            "--model-paths",
            str(checkpoint),
            "-o",
            str(output_path),
            "--accelerator",
            "cpu",
            "-q",
        ]
        if rdkit_2d:
            command.extend(["--molecule-featurizers", "rdkit_2d"])
        subprocess.run(command, check=True, capture_output=True, text=True)
        output = pd.read_csv(output_path)
    columns = [column for column in output if column.startswith(target_prefix)]
    if len(columns) != 1:
        raise KeyError(f"Expected one {target_prefix!r} output, found {columns}")
    return output[columns[0]].to_numpy(dtype=float)


def _tree_predictions(
    descriptor_frame: pd.DataFrame,
    artifact_path: Path,
    targets: list[str],
) -> dict[str, np.ndarray]:
    artifact: dict[str, Any] = joblib.load(artifact_path)
    x = descriptor_frame[artifact["descriptor_columns"]].to_numpy(dtype=np.float32)
    return {target: artifact["models"][target].predict(x) for target in targets}


def teacher_features(
    canonical_smiles: list[str], model_dir: str | Path = DEFAULT_MODEL_DIR
) -> tuple[pd.DataFrame, np.ndarray]:
    """Return deterministic structure descriptors and distilled-physics features."""
    model_dir = Path(model_dir)
    records = pd.DataFrame(
        {
            "molecule_id": [f"query_{index}" for index in range(len(canonical_smiles))],
            "canonical_smiles": canonical_smiles,
        }
    )
    descriptors = rdkit_descriptor_frame(records)
    qm = _chemprop_prediction(
        canonical_smiles,
        model_dir / "combisolv_qm.pt",
        "target_qm_water",
        rdkit_2d=True,
    )
    smd = _chemprop_prediction(
        canonical_smiles,
        model_dir / "molsolv_smd.pt",
        "target_smd_water",
        rdkit_2d=False,
    )
    abraham_targets = [f"abraham_{name}" for name in "esabl"]
    abraham = _tree_predictions(descriptors, model_dir / "abraham_teacher.joblib", abraham_targets)
    openff = _tree_predictions(
        descriptors,
        model_dir / "openff_teacher.joblib",
        ["openff23_dg", "openff23_exp_residual"],
    )
    implicit = _tree_predictions(
        descriptors,
        model_dir / "implicit_teacher.joblib",
        ["gbn2_alchemical_dg", "gbn2_exp_residual"],
    )
    confsolv_targets = [
        "confsolv_gas_conformer_correction",
        "confsolv_solution_conformer_correction",
        "confsolv_hydration_conformer_correction",
        "confsolv_water_gsolv_std",
        "confsolv_water_response_mean",
        "confsolv_water_response_std",
    ]
    confsolv = _tree_predictions(
        descriptors, model_dir / "confsolv_teacher.joblib", confsolv_targets
    )
    physical = np.column_stack(
        [
            qm,
            *(abraham[target] for target in abraham_targets),
            openff["openff23_dg"] + openff["openff23_exp_residual"],
            implicit["gbn2_alchemical_dg"] + implicit["gbn2_exp_residual"],
            smd,
            *(confsolv[target] for target in confsolv_targets),
        ]
    )
    return descriptors, physical


def predict_smiles(
    smiles: list[str], model_dir: str | Path = DEFAULT_MODEL_DIR
) -> tuple[np.ndarray, np.ndarray]:
    """Predict hydration free energy and seed-ensemble spread from SMILES only."""
    model_dir = Path(model_dir)
    artifact: dict[str, Any] = joblib.load(model_dir / "head.joblib")
    canonical = [canonicalize(value)[0] for value in smiles]
    descriptors, physical = teacher_features(canonical, model_dir)
    x = np.column_stack(
        [descriptors[artifact["descriptor_columns"]].to_numpy(dtype=float), physical]
    )
    if x.shape[1] != len(artifact["descriptor_columns"]) + len(artifact["teacher_columns"]):
        raise AssertionError("Inference feature schema does not match the fitted head")
    member_predictions = np.column_stack([model.predict(x) for model in artifact["models"]])
    return member_predictions.mean(axis=1), member_predictions.std(axis=1)
