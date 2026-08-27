"""Structure-to-response surrogates bundled with SolvAI.

These models predict physical response priors from molecular structure. They do
not run SMD, molecular dynamics, PIMD or any other simulation at inference.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd


def _chemprop_prediction(
    smiles: list[str], checkpoint: Path, target_prefix: str, *, rdkit_2d: bool
) -> np.ndarray:
    executable = shutil.which("chemprop")
    if executable is None:
        raise RuntimeError("chemprop executable not found; run `uv sync --extra inference`")
    with tempfile.TemporaryDirectory(prefix="solvai_") as directory:
        temporary = Path(directory)
        input_path = temporary / "input.csv"
        output_path = temporary / "output.csv"
        pd.DataFrame({"smiles": smiles}).to_csv(input_path, index=False)
        command = [
            executable,
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
    descriptors: pd.DataFrame, artifact_path: Path, targets: list[str]
) -> dict[str, np.ndarray]:
    artifact: dict[str, Any] = joblib.load(artifact_path)
    values = descriptors[artifact["descriptor_columns"]].to_numpy(dtype=np.float32)
    return {target: artifact["models"][target].predict(values) for target in targets}


def response_features(smiles: list[str], descriptors: pd.DataFrame, model_dir: Path) -> np.ndarray:
    """Predict the 15 frozen physical-response features from structure."""
    qm = _chemprop_prediction(
        smiles, model_dir / "combisolv_qm.pt", "target_qm_water", rdkit_2d=True
    )
    smd = _chemprop_prediction(
        smiles, model_dir / "molsolv_smd.pt", "target_smd_water", rdkit_2d=False
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
    return np.column_stack(
        [
            qm,
            *(abraham[target] for target in abraham_targets),
            openff["openff23_dg"] + openff["openff23_exp_residual"],
            implicit["gbn2_alchemical_dg"] + implicit["gbn2_exp_residual"],
            smd,
            *(confsolv[target] for target in confsolv_targets),
        ]
    )
