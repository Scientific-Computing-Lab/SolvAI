"""Build structure-linked ConfSolv water-response teacher targets.

ConfSolv provides conformer-resolved DFT/COSMO-RS solution thermodynamics.  We
collapse the water rows into a physically coherent solute-level hierarchy:
gas and solution conformational free-energy corrections, the resulting
Boltzmann-ensemble hydration free energy, and distributional summaries of the
conformer-specific COSMO-RS response.  All ARROW benchmark connectivities are
removed before the table is written.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from rdkit import Chem

R_KJ_MOL_K = 0.00831446261815324
TEMPERATURE_K = 298.0
RT_KJ_MOL = R_KJ_MOL_K * TEMPERATURE_K
KJ_PER_KCAL = 4.184


def boltzmann_correction(relative_energy: np.ndarray) -> float:
    values = np.asarray(relative_energy, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    minimum = float(values.min())
    return float(minimum - RT_KJ_MOL * np.log(np.exp(-(values - minimum) / RT_KJ_MOL).sum()))


def summarize_group(group: pd.DataFrame) -> dict[str, object]:
    gas = group["G(gas)"].to_numpy(dtype=np.float64)
    solution = group["G(solution)"].to_numpy(dtype=np.float64)
    solvation = group["G(solvation)"].to_numpy(dtype=np.float64)
    valid = np.isfinite(gas) & np.isfinite(solution) & np.isfinite(solvation)
    gas = gas[valid]
    solution = solution[valid]
    solvation = solvation[valid]
    if not len(gas):
        return {"confsolv_conformer_count": 0}

    gas_relative = gas - gas.min()
    solution_relative = solution - solution.min()
    gas_correction = boltzmann_correction(gas_relative)
    solution_correction = boltzmann_correction(solution_relative)
    ensemble_hydration = solution.min() + solution_correction - gas.min() - gas_correction
    gas_lec = int(np.argmin(gas))
    solution_lec = int(np.argmin(solution))
    response = solution_relative - gas_relative
    weights_gas = np.exp(-gas_relative / RT_KJ_MOL)
    weights_gas /= weights_gas.sum()
    weights_solution = np.exp(-solution_relative / RT_KJ_MOL)
    weights_solution /= weights_solution.sum()

    result = {
        "confsolv_conformer_count": len(gas),
        "confsolv_water_dg_ensemble": float(ensemble_hydration),
        "confsolv_water_gsolv_gas_lec": float(solvation[gas_lec]),
        "confsolv_water_gsolv_solution_lec": float(solvation[solution_lec]),
        "confsolv_water_gsolv_gas_weighted": float(np.dot(weights_gas, solvation)),
        "confsolv_water_gsolv_solution_weighted": float(np.dot(weights_solution, solvation)),
        "confsolv_gas_conformer_correction": float(gas_correction),
        "confsolv_solution_conformer_correction": float(solution_correction),
        "confsolv_hydration_conformer_correction": float(ensemble_hydration - solvation[gas_lec]),
        "confsolv_water_gsolv_mean": float(solvation.mean()),
        "confsolv_water_gsolv_std": float(solvation.std()),
        "confsolv_water_gsolv_min": float(solvation.min()),
        "confsolv_water_gsolv_max": float(solvation.max()),
        "confsolv_water_gsolv_q10": float(np.quantile(solvation, 0.10)),
        "confsolv_water_gsolv_q50": float(np.quantile(solvation, 0.50)),
        "confsolv_water_gsolv_q90": float(np.quantile(solvation, 0.90)),
        "confsolv_water_response_mean": float(response.mean()),
        "confsolv_water_response_std": float(response.std()),
        "confsolv_water_response_min": float(response.min()),
        "confsolv_water_response_max": float(response.max()),
        "confsolv_water_response_q10": float(np.quantile(response, 0.10)),
        "confsolv_water_response_q50": float(np.quantile(response, 0.50)),
        "confsolv_water_response_q90": float(np.quantile(response, 0.90)),
    }
    for name in tuple(result):
        if name != "confsolv_conformer_count":
            result[name] = float(result[name]) / KJ_PER_KCAL
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--water",
        type=Path,
        default=ROOT / "data/processed/confsolv_water_raw.parquet",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=ROOT / "data/processed/confsolv_structure_mapping.parquet",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/processed/confsolv_water_nonbenchmark.parquet",
    )
    args = parser.parse_args()

    water = pd.read_parquet(args.water)
    required = {"mol_id", "conf_id", "solvent", "G(gas)", "G(solution)", "G(solvation)"}
    missing = required - set(water)
    if missing:
        raise KeyError(f"Missing ConfSolv columns: {sorted(missing)}")
    if not water.solvent.eq("H2O").all():
        raise ValueError("Input contains non-water ConfSolv rows")

    rows = []
    grouped = water.groupby("mol_id", sort=True, observed=True)
    for index, (mol_id, group) in enumerate(grouped, start=1):
        rows.append({"confsolv_id": str(mol_id), **summarize_group(group)})
        if index % 5000 == 0:
            print(f"Aggregated {index}/{grouped.ngroups} solutes", flush=True)
    physics = pd.DataFrame(rows)
    mapping = pd.read_parquet(args.mapping)
    source = physics.merge(mapping, on="confsolv_id", how="left", validate="one_to_one")
    source = source[source.mapping_status.eq("resolved")].copy()
    source = source.dropna(subset=["canonical_smiles", "connectivity_key"])

    neutral = []
    heavy_atoms = []
    for smiles in source.canonical_smiles:
        mol = Chem.MolFromSmiles(smiles)
        neutral.append(mol is not None and Chem.GetFormalCharge(mol) == 0)
        heavy_atoms.append(np.nan if mol is None else mol.GetNumHeavyAtoms())
    source["formal_charge_zero"] = neutral
    source["heavy_atom_count"] = heavy_atoms
    source = source[source.formal_charge_zero].copy()

    benchmark = pd.read_parquet(ROOT / "data/processed/arrow_solvation_master.parquet")
    benchmark_keys = set(
        benchmark.loc[benchmark.solvent.eq("water"), "inchi_connectivity_key"].astype(str)
    )
    overlaps = set(source.connectivity_key.astype(str)) & benchmark_keys
    source = source[~source.connectivity_key.astype(str).isin(benchmark_keys)].copy()

    numeric = [
        column
        for column in source
        if column.startswith("confsolv_") and pd.api.types.is_numeric_dtype(source[column])
    ]
    aggregation: dict[str, str] = {
        "confsolv_id": "first",
        "canonical_smiles": "first",
        "inchi_key": "first",
        "heavy_atom_count": "first",
        **{column: "median" for column in numeric},
    }
    source = (
        source.groupby("connectivity_key", as_index=False, sort=True)
        .agg(aggregation)
        .reset_index(drop=True)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    source.to_parquet(args.output, index=False)
    metadata = {
        "source": "ConfSolv DOI 10.5281/zenodo.8292520",
        "method": "TPSS-D3(BJ)/def2-TZVP + GFN2-xTB RRHO + COSMO-RS BP-TZVP",
        "temperature_k": TEMPERATURE_K,
        "archive_energy_units": "kJ/mol",
        "output_energy_units": "kcal/mol",
        "raw_water_rows": len(water),
        "raw_water_solutes": int(water.mol_id.nunique()),
        "resolved_neutral_unique_connectivities_after_exclusion": len(source),
        "benchmark_connectivity_overlaps_removed": len(overlaps),
        "target_columns": numeric,
        "inference_simulation": False,
        "role": "training-time privileged water-response supervision",
    }
    metadata_path = args.output.with_name("confsolv_water_teacher_metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
