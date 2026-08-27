"""Deterministic molecular identities and two-dimensional descriptors."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator


def canonicalize(smiles: str) -> tuple[str, str, str]:
    """Return canonical isomeric SMILES, full InChIKey and connectivity block."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    key = Chem.MolToInchiKey(molecule)
    return canonical, key, key.split("-")[0]


def descriptor_frame(smiles: list[str]) -> pd.DataFrame:
    """Compute the descriptor schema used by the frozen SolvAI artifact."""
    fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rows: list[dict[str, Any]] = []
    for value in smiles:
        canonical, _, _ = canonicalize(value)
        molecule = Chem.MolFromSmiles(canonical)
        if molecule is None:  # pragma: no cover - guarded by canonicalize
            raise ValueError(f"RDKit could not parse canonical SMILES: {canonical}")
        row: dict[str, Any] = {"canonical_smiles": canonical}
        for name, function in Descriptors._descList:
            try:
                result = float(function(molecule))
                row[f"rdkit__{name}"] = result if math.isfinite(result) else np.nan
            except (ValueError, OverflowError, RuntimeError, ZeroDivisionError):
                row[f"rdkit__{name}"] = np.nan
        bits = fingerprint.GetFingerprintAsNumPy(molecule)
        row.update({f"morgan2__{index:04d}": int(bit) for index, bit in enumerate(bits)})
        rows.append(row)
    return pd.DataFrame(rows)
