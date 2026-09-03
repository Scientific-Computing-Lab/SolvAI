#!/usr/bin/env python3
"""Build the frozen Gate-1 time-series and molecule-feature tables."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem

from solv_ai.features import descriptor_frame
from solv_ai.teachers import response_features

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai_v3/data/derived"
PHASE2 = ROOT / "active_solvai/results/phase2"
MODEL = ROOT / "models/final"

STRUCTURE_COLUMNS = {
    "structure_molecular_weight": "rdkit__MolWt",
    "structure_heavy_atom_count": "rdkit__HeavyAtomCount",
    "structure_tpsa": "rdkit__TPSA",
    "structure_mol_logp": "rdkit__MolLogP",
    "structure_hbond_donors": "rdkit__NumHDonors",
    "structure_hbond_acceptors": "rdkit__NumHAcceptors",
    "structure_rotatable_bonds": "rdkit__NumRotatableBonds",
    "structure_ring_count": "rdkit__RingCount",
    "structure_fraction_csp3": "rdkit__FractionCSP3",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_energy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [str(column).strip() for column in frame.columns]
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    required = ["Time", "Temp", "Density", "dHdL"]
    missing = [column for column in required if column not in numeric]
    if missing:
        raise AssertionError(f"{path} lacks {missing}")
    if len(numeric) != 51:
        raise AssertionError(f"{path} has {len(numeric)} frames, expected 51")
    post = numeric[numeric["Time"] > 1e-9].copy().reset_index(drop=True)
    if len(post) != 50 or not np.isfinite(post[required]).all().all():
        raise AssertionError(f"{path} does not contain 50 finite post-initial records")
    post["frame_index"] = np.arange(1, 51)
    post["block_0p5ps"] = (post["frame_index"] - 1) // 5 + 1
    post["block_1ps"] = (post["frame_index"] - 1) // 10 + 1
    return post


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_paths = [
        PHASE2 / "dense_responses_calibration.parquet",
        PHASE2 / "dense_responses_prospective.parquet",
    ]
    responses = pd.concat([pd.read_parquet(path) for path in source_paths], ignore_index=True)
    responses = responses.rename(columns={"lambda": "lambda_value"})
    if len(responses) != 180 or responses.groupby("molecule_id").size().ne(15).any():
        raise AssertionError("Gate-1 source is not a complete 12 x 15 grid")

    time_rows: list[pd.DataFrame] = []
    raw_hashes: dict[str, str] = {}
    for row in responses.itertuples(index=False):
        path = Path(str(row.energy_file))
        observed_hash = sha256(path)
        if observed_hash != row.energy_sha256:
            raise AssertionError(f"energy hash mismatch: {path}")
        raw_hashes[str(path)] = observed_hash
        frame = read_energy(path)
        frame = frame.rename(
            columns={
                "Time": "time_ps",
                "Temp": "temperature_k",
                "Density": "density_g_cm3",
                "dHdL": "dhdl_kcal_mol",
            }
        )
        frame.insert(0, "lambda", float(row.lambda_value))
        frame.insert(0, "canonical_smiles", row.canonical_smiles)
        frame.insert(0, "molecule_name", row.molecule_name)
        frame.insert(0, "molecule_id", row.molecule_id)
        frame["energy_file"] = str(path)
        frame["energy_sha256"] = observed_hash
        time_rows.append(
            frame[
                [
                    "molecule_id",
                    "molecule_name",
                    "canonical_smiles",
                    "lambda",
                    "frame_index",
                    "block_0p5ps",
                    "block_1ps",
                    "time_ps",
                    "dhdl_kcal_mol",
                    "temperature_k",
                    "density_g_cm3",
                    "energy_file",
                    "energy_sha256",
                ]
            ]
        )
    time_series = pd.concat(time_rows, ignore_index=True).sort_values(
        ["molecule_id", "lambda", "frame_index"]
    )
    time_path = OUT / "gate1_time_series.parquet"
    time_series.to_parquet(time_path, index=False)

    molecules = (
        responses[["molecule_id", "molecule_name", "canonical_smiles", "functional_group_family"]]
        .drop_duplicates()
        .sort_values("molecule_id")
        .reset_index(drop=True)
    )
    smiles = molecules.canonical_smiles.tolist()
    descriptors = descriptor_frame(smiles)
    physical = response_features(smiles, descriptors, MODEL)
    head = joblib.load(MODEL / "head.joblib")
    response_columns = [f"response__{name}" for name in head["teacher_columns"]]
    if physical.shape != (len(molecules), 15) or len(response_columns) != 15:
        raise AssertionError("frozen SolvAI response schema is not 15-dimensional")
    features = molecules.copy()
    for output_name, descriptor_name in STRUCTURE_COLUMNS.items():
        features[output_name] = descriptors[descriptor_name].to_numpy(float)
    features["structure_formal_charge"] = [
        int(Chem.GetFormalCharge(Chem.MolFromSmiles(value))) for value in smiles
    ]
    for index, name in enumerate(response_columns):
        features[name] = physical[:, index]
    feature_path = OUT / "gate1_molecule_features.parquet"
    features.to_parquet(feature_path, index=False)

    manifest = {
        "schema_version": 1,
        "molecules": int(features.molecule_id.nunique()),
        "windows": int(responses.shape[0]),
        "post_initial_frames": int(len(time_series)),
        "lambda_values": sorted(responses["lambda_value"].unique().tolist()),
        "structure_columns": [*STRUCTURE_COLUMNS, "structure_formal_charge"],
        "response_columns": response_columns,
        "source_tables": {str(path): sha256(path) for path in source_paths},
        "raw_energy_files": raw_hashes,
        "outputs": {
            str(time_path): sha256(time_path),
            str(feature_path): sha256(feature_path),
        },
        "independence_warning": (
            "One 5-ps stream per molecule-window; blocks are development diagnostics, not replicas."
        ),
    }
    manifest_path = OUT / "gate1_dataset_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in manifest.items() if key != "raw_energy_files"}, indent=2))


if __name__ == "__main__":
    main()
