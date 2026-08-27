"""Inventory and normalize external solvation/physics supervision sources."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd
from arrow_distill.data import ROOT
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def normalize_smiles(smiles: object) -> tuple[str | None, str | None, str | None, int | None]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None, None, None, None
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, None, None, None
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    key = Chem.MolToInchiKey(molecule)
    return canonical, key, key[:14], Chem.GetFormalCharge(molecule)


def map_smiles(
    values: pd.Series,
) -> dict[str, tuple[str | None, str | None, str | None, int | None]]:
    return {value: normalize_smiles(value) for value in values.dropna().astype(str).unique()}


def add_identity(frame: pd.DataFrame, column: str, prefix: str) -> pd.DataFrame:
    mapping = map_smiles(frame[column])
    identity = frame[column].map(mapping)
    identity = identity.map(
        lambda item: item if isinstance(item, tuple) else (None, None, None, None)
    )
    frame[f"{prefix}_canonical_smiles"] = [item[0] for item in identity]
    frame[f"{prefix}_inchi_key"] = [item[1] for item in identity]
    frame[f"{prefix}_connectivity_key"] = [item[2] for item in identity]
    frame[f"{prefix}_formal_charge"] = [item[3] for item in identity]
    return frame


def git_commit(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except subprocess.CalledProcessError:
        return None


def value_count_solv_tum(path: Path) -> tuple[int, int]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecules = 0
    pairs = 0
    for molecule in supplier:
        if molecule is None:
            continue
        molecules += 1
        pairs += sum(name.startswith("logK (") for name in molecule.GetPropNames())
    return molecules, pairs


def benchmark_overlap(frame: pd.DataFrame, key: str, benchmark_keys: set[str]) -> int:
    return int(frame[key].isin(benchmark_keys).sum())


def main() -> None:
    processed = ROOT / "data/processed"
    processed.mkdir(parents=True, exist_ok=True)
    catalog_dir = ROOT / "data/catalog"
    catalog_dir.mkdir(parents=True, exist_ok=True)

    master = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = master[master.solvent.eq("water")].copy()
    benchmark_keys = set(benchmark.inchi_connectivity_key)
    des_water = pd.read_parquet(processed / "des370k_water_response_nonbenchmark.parquet")
    des_metadata = json.loads((processed / "des370k_water_teacher_metadata.json").read_text())
    sampl4_physics = pd.read_parquet(processed / "sampl4_physics_nonbenchmark.parquet")
    sampl4_metadata = json.loads((processed / "sampl4_physics_teacher_metadata.json").read_text())
    bannan_curves = pd.read_parquet(processed / "bannan_alchemical_curve_nonbenchmark.parquet")
    bannan_metadata = json.loads(
        (processed / "bannan_alchemical_curve_teacher_metadata.json").read_text()
    )
    reft = pd.read_csv(ROOT / "data/external/zenodo/20699238/hydration_free_enenergies.csv")
    openfe_diagnostics = pd.read_parquet(processed / "openfe_diagnostics_nonbenchmark.parquet")
    openfe_metadata = json.loads(
        (processed / "openfe_diagnostics_teacher_metadata.json").read_text()
    )
    nqeliq = pd.read_parquet(processed / "nqeliq_nonbenchmark.parquet")
    nqeliq_metadata = json.loads((processed / "nqeliq_teacher_metadata.json").read_text())
    molsolv_smd = pd.read_parquet(processed / "molsolv_smd_water_nonbenchmark.parquet")
    molsolv_metadata = json.loads((processed / "molsolv_smd_teacher_metadata.json").read_text())
    phase_space = pd.read_parquet(processed / "phase_space_dynamic_nonbenchmark.parquet")
    phase_space_metadata = json.loads(
        (processed / "phase_space_dynamic_teacher_metadata.json").read_text()
    )
    gnequip = pd.read_parquet(processed / "gnequip_solvation_nonbenchmark.parquet")
    gnequip_metadata = json.loads(
        (processed / "gnequip_solvation_teacher_metadata.json").read_text()
    )
    mlff_hfe = pd.read_parquet(processed / "mlff_hfe_nonbenchmark.parquet")
    mlff_hfe_metadata = json.loads((processed / "mlff_hfe_teacher_metadata.json").read_text())
    confsolv_water = pd.read_parquet(processed / "confsolv_water_nonbenchmark.parquet")
    confsolv_metadata = json.loads((processed / "confsolv_water_teacher_metadata.json").read_text())

    combi_qm_path = ROOT / "data/external/combisolv_qm_original.txt"
    combi_qm_cache = processed / "combisolv_qm_all.parquet"
    if combi_qm_cache.exists():
        combi_qm = pd.read_parquet(combi_qm_cache)
    else:
        combi_qm = pd.read_csv(combi_qm_path)
        combi_qm = add_identity(combi_qm, "mol solute", "solute")
        combi_qm = add_identity(combi_qm, "mol solvent", "solvent")
        combi_qm = combi_qm.rename(columns={"target Gsolv kcal": "delta_g_solv_qm"})
        combi_qm["benchmark_solute_overlap"] = combi_qm.solute_connectivity_key.isin(benchmark_keys)
        combi_qm["source_method"] = "COSMOtherm (CombiSolv-QM)"
        combi_qm["source_doi"] = "10.1016/j.cej.2021.129307"
        combi_qm.to_parquet(combi_qm_cache, index=False)
    strict_qm = combi_qm[
        combi_qm.solute_canonical_smiles.notna()
        & combi_qm.solvent_canonical_smiles.notna()
        & ~combi_qm.benchmark_solute_overlap
    ].copy()
    strict_qm.to_parquet(processed / "combisolv_qm_nonbenchmark.parquet", index=False)
    water_qm = strict_qm[strict_qm.solvent_canonical_smiles.eq("O")].copy()
    water_qm.to_parquet(processed / "combisolv_qm_water_nonbenchmark.parquet", index=False)

    combi_exp_path = ROOT / "data/external/combisolv_exp_original.xlsx"
    combi_exp = pd.read_excel(combi_exp_path, sheet_name="data")
    combi_exp = add_identity(combi_exp, "smiles_solute", "solute")
    combi_exp = add_identity(combi_exp, "smiles_solvent", "solvent")
    combi_exp["benchmark_solute_overlap"] = combi_exp.solute_connectivity_key.isin(benchmark_keys)
    strict_exp = combi_exp[
        combi_exp.solute_canonical_smiles.notna()
        & combi_exp.solvent_canonical_smiles.notna()
        & ~combi_exp.benchmark_solute_overlap
    ].copy()
    strict_exp.to_parquet(processed / "combisolv_exp_all_nonbenchmark.parquet", index=False)

    g4_path = ROOT / "data/external/g4mp2_solvation/train.parquet"
    g4_cache = processed / "g4mp2_solvation_all.parquet"
    if g4_cache.exists():
        g4 = pd.read_parquet(g4_cache)
    else:
        g4 = pd.read_parquet(g4_path)
        g4 = add_identity(g4, "smiles_0", "solute")
        g4["benchmark_solute_overlap"] = g4.solute_connectivity_key.isin(benchmark_keys)
        g4["source_method"] = "DFT continuum solvation on G4MP2/QM9 geometries"
        g4["source_doi"] = "10.18126/jos5-wj65"
        g4.to_parquet(g4_cache, index=False)
    strict_g4 = g4[
        g4.solute_canonical_smiles.notna() & g4.sol_water.notna() & ~g4.benchmark_solute_overlap
    ].copy()
    strict_g4.to_parquet(processed / "g4mp2_solvation_nonbenchmark.parquet", index=False)

    guthrie_path = (
        ROOT / "data/external/repos/jazzy/data/free_energy_guthrie/02-guthrie_curated.csv"
    )
    guthrie = pd.read_csv(guthrie_path)
    guthrie = add_identity(guthrie, "mol", "solute")
    guthrie["benchmark_solute_overlap"] = guthrie.solute_connectivity_key.isin(benchmark_keys)
    guthrie.to_parquet(processed / "guthrie_curated_identity.parquet", index=False)

    sampl_path = ROOT / "data/external/repos/SAMPL1/SAMPL1_Dataset.csv"
    sampl = pd.read_csv(sampl_path)
    sampl = add_identity(sampl, "SMILES", "solute")
    sampl["benchmark_solute_overlap"] = sampl.solute_connectivity_key.isin(benchmark_keys)
    strict_sampl = sampl[
        sampl.solute_canonical_smiles.notna() & ~sampl.benchmark_solute_overlap
    ].copy()
    strict_sampl.to_parquet(processed / "sampl1_nonbenchmark.parquet", index=False)

    solvatum_path = ROOT / "data/external/repos/solvatum/solvatum/data/solvatum.sdf"
    solvatum_molecules, solvatum_pairs = value_count_solv_tum(solvatum_path)
    thermo_root = ROOT / "data/external/repos/thermo-length-learning"
    thermo_molecules = sum(
        1
        for line in (thermo_root / "database.txt").read_text().splitlines()
        if line and not line.startswith("#")
    )
    thermo_pairs = len(list((thermo_root / "results").glob("lig*to*")))
    probe_path = ROOT / "results/probe_results.parquet"
    probe_count = 0
    if probe_path.exists():
        probe_count = int(pd.read_parquet(probe_path).success.astype("boolean").fillna(False).sum())

    catalog = [
        {
            "source": "Freecurve ARROW hydration benchmark",
            "molecules_or_samples": "85 molecules / 125 molecule-solvent rows",
            "targets": "experiment; classical ARROW; PIMD4; PIMD8; NQE and experimental residuals",
            "physics_content": "high-fidelity alchemical free energies and bead convergence",
            "water_specific": "Primary 85 water; 40 cyclohexane auxiliary",
            "data_class": "experimental + classical MD + PIMD",
            "downloaded": True,
            "usable_for_training": "Yes; outer-fold masked labels only",
            "canonical_url": "https://doi.org/10.1038/s41467-022-28041-0",
            "doi": "10.1038/s41467-022-28041-0",
            "license": "Article supplementary data; check publisher terms",
            "exact_data_fields": ";".join(master.columns),
            "benchmark_overlap": 85,
            "pretraining_use": "Primary downstream benchmark",
            "privileged_supervision": "Classical/PIMD/NQE/residual heads",
            "leakage_risk": "Experimental target must be outer-fold isolated",
            "local_path": "data/processed/arrow_solvation_master.parquet",
        },
        {
            "source": "Freecurve short PIMD2 lambda=0.5 probes",
            "molecules_or_samples": f"{probe_count} completed of 76 configured",
            "targets": "dH/dlambda and component/time-series/shell/bead statistics",
            "physics_content": "short explicit-water PIMD2 response fingerprints",
            "water_specific": "Yes",
            "data_class": "PIMD2 training-only teacher",
            "downloaded": True,
            "usable_for_training": "Yes; never an inference input",
            "canonical_url": "local generated data",
            "doi": None,
            "license": "Internal Freecurve",
            "exact_data_fields": "component energy and dH/dlambda summaries; solvent-shell geometry",
            "benchmark_overlap": probe_count,
            "pretraining_use": "Physics reconstruction heads",
            "privileged_supervision": "Yes",
            "leakage_risk": "Safe only as training labels inside each outer fold",
            "local_path": "results/probe_results.parquet",
        },
        {
            "source": "FreeSolv v0.52",
            "molecules_or_samples": "642 molecules",
            "targets": "experimental and GAFF hydration free energy with uncertainty",
            "physics_content": "experimental hydration plus classical alchemical estimates",
            "water_specific": "Yes",
            "data_class": "experimental + MD",
            "downloaded": True,
            "usable_for_training": "Yes after benchmark exclusion",
            "canonical_url": "https://github.com/MobleyLab/FreeSolv",
            "doi": "10.1021/ci3001277; 10.1021/acs.jced.7b00104",
            "license": "CC-BY-4.0",
            "exact_data_fields": "SMILES; experimental/calculated hydration free energy; uncertainties; references",
            "benchmark_overlap": 80,
            "pretraining_use": "Experimental hydration and classical teacher",
            "privileged_supervision": "Calculated hydration value",
            "leakage_risk": "80 exact benchmark connectivities globally excluded",
            "local_path": "data/raw/public/freesolv_database.json",
        },
        {
            "source": "CombiSolv-Exp-8780",
            "molecules_or_samples": f"{len(combi_exp):,} solute-solvent records; {combi_exp.solute_connectivity_key.nunique():,} solutes",
            "targets": "experimental solvation free energy and replicate spread",
            "physics_content": "solvent-conditioned experimental response",
            "water_specific": f"No; {int(combi_exp.solvent_canonical_smiles.eq('O').sum()):,} water rows",
            "data_class": "experimental",
            "downloaded": True,
            "usable_for_training": "Yes after solute-level benchmark exclusion",
            "canonical_url": "https://ars.els-cdn.com/content/image/1-s2.0-S1385894721008925-mmc2.xlsx",
            "doi": "10.1016/j.cej.2021.129307",
            "license": "Supplementary-file reuse terms; no standalone license",
            "exact_data_fields": ";".join(combi_exp.columns[:12]),
            "benchmark_overlap": benchmark_overlap(
                combi_exp, "solute_connectivity_key", benchmark_keys
            ),
            "pretraining_use": "Solvent-conditioned representation / water target",
            "privileged_supervision": "No",
            "leakage_risk": "Water subset overlaps benchmark; all benchmark solutes removed",
            "local_path": "data/processed/combisolv_exp_all_nonbenchmark.parquet",
        },
        {
            "source": "CombiSolv-QM",
            "molecules_or_samples": f"{len(combi_qm):,} pairs; {combi_qm.solute_connectivity_key.nunique():,} solutes; {combi_qm.solvent_connectivity_key.nunique():,} solvents",
            "targets": "COSMOtherm solvation free energy",
            "physics_content": "quantum-continuum solvent response across 284 solvents",
            "water_specific": f"No; {int(combi_qm.solvent_canonical_smiles.eq('O').sum()):,} water rows",
            "data_class": "QM/continuum",
            "downloaded": True,
            "usable_for_training": "Yes after strict solute exclusion",
            "canonical_url": "https://ars.els-cdn.com/content/image/1-s2.0-S1385894721008925-mmc1.txt",
            "doi": "10.1016/j.cej.2021.129307",
            "license": "Supplementary-file reuse terms; no standalone license",
            "exact_data_fields": "solvent SMILES; solute SMILES; Gsolv (kcal/mol)",
            "benchmark_overlap": benchmark_overlap(
                combi_qm, "solute_connectivity_key", benchmark_keys
            ),
            "pretraining_use": "Water-only teacher and solvent-conditioned pretraining",
            "privileged_supervision": "Yes: approximate physics target",
            "leakage_risk": "Benchmark solutes removed from every solvent row in strict files",
            "local_path": "data/processed/combisolv_qm_nonbenchmark.parquet",
        },
        {
            "source": "G4MP2 solvation/QM9 (Foundry-ML)",
            "molecules_or_samples": f"{len(g4):,} molecules",
            "targets": "five solvent energies; charges; dipole; polarizability; orbital and thermochemical properties",
            "physics_content": "DFT continuum solvation plus quantum molecular observables",
            "water_specific": "Includes water and four organic solvents",
            "data_class": "QM/DFT",
            "downloaded": True,
            "usable_for_training": "Yes after strict benchmark exclusion",
            "canonical_url": "https://huggingface.co/datasets/foundry-ml/foundry_g4mp2_solvation_v1-2",
            "doi": "10.18126/jos5-wj65",
            "license": "CC-BY-4.0",
            "exact_data_fields": ";".join(g4.columns[:36]),
            "benchmark_overlap": benchmark_overlap(g4, "solute_connectivity_key", benchmark_keys),
            "pretraining_use": "Multi-head quantum and solvation encoder",
            "privileged_supervision": "Yes: solvation/charge/dipole/polarizability heads",
            "leakage_risk": "All exact benchmark connectivities removed in strict file",
            "local_path": "data/processed/g4mp2_solvation_nonbenchmark.parquet",
        },
        {
            "source": "GuthrieSolv",
            "molecules_or_samples": "53,895 raw literature records; 3,316 Jazzy-curated rows",
            "targets": "heterogeneous hydration-related measurements",
            "physics_content": "experimental source records and uncertainty/context",
            "water_specific": "Mostly hydration; raw table mixes properties/units",
            "data_class": "experimental",
            "downloaded": True,
            "usable_for_training": "Curated subsets only; not used blindly",
            "canonical_url": "https://github.com/MobleyLab/GuthrieSolv",
            "doi": "10.5281/zenodo.1101258",
            "license": "CC-BY-4.0",
            "exact_data_fields": "SMILES; name; process; value; uncertainty; unit; temperature; reference",
            "benchmark_overlap": benchmark_overlap(
                guthrie, "solute_connectivity_key", benchmark_keys
            ),
            "pretraining_use": "Potential after process/unit curation",
            "privileged_supervision": "No",
            "leakage_risk": "Substantial FreeSolv/benchmark ancestry",
            "local_path": "data/external/repos/GuthrieSolv/guthrie_database.csv",
        },
        {
            "source": "SAMPL1 hydration challenge",
            "molecules_or_samples": f"{len(sampl)} molecules",
            "targets": "experimental hydration free energy and uncertainty",
            "physics_content": "blind-challenge experimental reference values and 3D structures",
            "water_specific": "Yes",
            "data_class": "experimental benchmark",
            "downloaded": True,
            "usable_for_training": "Yes after benchmark exclusion; small and out-of-domain",
            "canonical_url": "https://github.com/leelasd/SAMPL1",
            "doi": "10.1021/jp806724u",
            "license": "No repository license found",
            "exact_data_fields": "challenge key; name; Ghyd; uncertainty; SMILES; PDB",
            "benchmark_overlap": benchmark_overlap(
                sampl, "solute_connectivity_key", benchmark_keys
            ),
            "pretraining_use": "Hydration auxiliary target",
            "privileged_supervision": "No",
            "leakage_risk": "Exact benchmark connectivities removed",
            "local_path": "data/processed/sampl1_nonbenchmark.parquet",
        },
        {
            "source": "Solv@TUM",
            "molecules_or_samples": f"{solvatum_molecules:,} solutes / {solvatum_pairs:,} logK pairs",
            "targets": "gas-solvent partition coefficient converted to solvation free energy",
            "physics_content": "experimental multi-solvent response; dipole/polarizability",
            "water_specific": "No; non-aqueous",
            "data_class": "experimental",
            "downloaded": True,
            "usable_for_training": "Yes as solvent-conditioned auxiliary task",
            "canonical_url": "https://github.com/hille721/solvatum",
            "doi": "10.14459/2018mp1452571.001",
            "license": "CC-BY-SA-4.0",
            "exact_data_fields": "SMILES; InChI; logK by solvent; polarizability; dipole",
            "benchmark_overlap": "Not yet materialized at pair level",
            "pretraining_use": "Solvent-conditioned representation",
            "privileged_supervision": "Dipole/polarizability",
            "leakage_risk": "Exclude benchmark solutes before supervised use",
            "local_path": "data/external/repos/solvatum/solvatum/data/solvatum.sdf",
        },
        {
            "source": "Thermodynamic-length learning dataset",
            "molecules_or_samples": f"{thermo_molecules} FreeSolv molecules / {thermo_pairs} relative transformations",
            "targets": "relative hydration free energies and variances",
            "physics_content": "alchemical network variance and trajectory-derived features",
            "water_specific": "Yes",
            "data_class": "classical MD",
            "downloaded": True,
            "usable_for_training": "Potentially, after pair/trajectory parsing",
            "canonical_url": "https://github.com/choderalab/thermo-length-learning",
            "doi": "10.48550/arXiv.1906.08599",
            "license": "No repository license found",
            "exact_data_fields": "FreeSolv records; pair transformations; result pickles; variance analyses",
            "benchmark_overlap": "FreeSolv-derived; expected high",
            "pretraining_use": "Alchemical response/uncertainty auxiliary task",
            "privileged_supervision": "Yes",
            "leakage_risk": "Experimental labels excluded; physics features need identity audit",
            "local_path": "data/external/repos/thermo-length-learning",
        },
        {
            "source": "Jazzy",
            "molecules_or_samples": "292 Gerber fit molecules; 3,316 Guthrie validation rows",
            "targets": "hydration free energy and H-bond donor/acceptor strengths",
            "physics_content": "EEQ charges, atomic polarizabilities, polar/apolar decomposition",
            "water_specific": "Yes",
            "data_class": "deterministic QSPR",
            "downloaded": True,
            "usable_for_training": "Only label-free atomic strengths are clean headline features",
            "canonical_url": "https://github.com/AstraZeneca/jazzy",
            "doi": "10.1038/s41598-023-30089-x",
            "license": "Apache-2.0",
            "exact_data_fields": "atomic EEQ charge/polarizability; sdc/sdx/sa; fitted dGa/dGp/dGi",
            "benchmark_overlap": "Gerber/Guthrie overlap audited separately",
            "pretraining_use": "Cheap deterministic physics descriptors",
            "privileged_supervision": "No trajectory; physics-informed representation",
            "leakage_risk": "Published dG coefficients fit experimental data; exclude dG from headline",
            "local_path": "data/external/repos/jazzy",
        },
        {
            "source": "SolvBERT",
            "molecules_or_samples": "1,000,000 QM + 8,780 experimental pairs",
            "targets": "solvation free energy; solubility",
            "physics_content": "solvent-solute language-model pretraining",
            "water_specific": "No",
            "data_class": "pretrained architecture/data",
            "downloaded": True,
            "usable_for_training": "Data yes; no downloadable checkpoint found",
            "canonical_url": "https://github.com/su-group/SolvBERT",
            "doi": "10.26434/chemrxiv-2022-0hl5p",
            "license": "No repository license found",
            "exact_data_fields": "paired SMILES; Gsolv",
            "benchmark_overlap": "See CombiSolv rows",
            "pretraining_use": "Architecture/reference; data used directly",
            "privileged_supervision": "QM solvent response",
            "leakage_risk": "Original experimental split is not benchmark-safe",
            "local_path": "data/external/repos/SolvBERT",
        },
        {
            "source": "3DMRL",
            "molecules_or_samples": "Uses CombiSolv pretraining; checkpoint absent",
            "targets": "molecular interaction properties",
            "physics_content": "virtual 3D solute-solvent interaction environment",
            "water_specific": "No",
            "data_class": "3D pretrained architecture",
            "downloaded": True,
            "usable_for_training": "Architecture feasible; full retraining not first-line",
            "canonical_url": "https://github.com/Namkyeong/3DMRL",
            "doi": "10.48550/arXiv.2412.02957",
            "license": "No repository license found",
            "exact_data_fields": "code only; external CombiSolv files required",
            "benchmark_overlap": "Depends on CombiSolv preprocessing",
            "pretraining_use": "Candidate 3D interaction encoder",
            "privileged_supervision": "Virtual solvent geometries",
            "leakage_risk": "Must rebuild splits; published splits are unsafe here",
            "local_path": "data/external/repos/3DMRL",
        },
        {
            "source": "GeoMAW-Solv / Solvaformer",
            "molecules_or_samples": "Code; external BigSolDBv2/CombiSolv",
            "targets": "solubility and solvation tasks",
            "physics_content": "3D Equiformer/MPNN; optional learned charges",
            "water_specific": "No",
            "data_class": "3D/graph architecture",
            "downloaded": True,
            "usable_for_training": "Code usable; advertised checkpoints require W&B access",
            "canonical_url": "https://github.com/Sanofi-Public/geomaw-solv",
            "doi": None,
            "license": "Repository license file; component licenses vary",
            "exact_data_fields": "code and configs; no main pretrained checkpoint in repository",
            "benchmark_overlap": "Depends on external data",
            "pretraining_use": "Candidate 3D encoder",
            "privileged_supervision": "Geometry/charge auxiliaries",
            "leakage_risk": "Must rebuild splits",
            "local_path": "data/external/repos/geomaw-solv",
        },
        {
            "source": "DES370K / DES5M",
            "molecules_or_samples": (
                f"{des_metadata['water_dimer_geometries']:,} neutral water-dimer geometries / "
                f"{des_metadata['external_teacher_solutes']:,} strict solutes; "
                "370,000 total DES370K geometries"
            ),
            "targets": "CCSD(T), SNS-MP2, MP2, HF, SAPT0 interaction energies",
            "physics_content": (
                "water-interaction distributions plus SAPT electrostatic, exchange, "
                "induction, and dispersion components"
            ),
            "water_specific": "Yes for the extracted teacher subset",
            "data_class": "QM dimers",
            "downloaded": True,
            "usable_for_training": "Yes; aggregated structure-only teacher built",
            "canonical_url": "https://zenodo.org/records/5676266; https://zenodo.org/records/5706002",
            "doi": "10.1038/s41597-021-00833-x",
            "license": "DESRES Data Sets License Agreement",
            "exact_data_fields": (
                "solute/water SMILES; sampled dimer geometries; CBS CCSD(T); "
                "SAPT electrostatic/exchange/induction/dispersion"
            ),
            "benchmark_overlap": (
                f"{des_metadata['benchmark_connectivity_overlap_removed']} connectivities removed; "
                "0 in strict teacher source"
            ),
            "pretraining_use": "Water interaction-response auxiliary heads",
            "privileged_supervision": "Yes",
            "leakage_risk": "All exact benchmark connectivities globally excluded",
            "local_path": "data/processed/des370k_water_response_nonbenchmark.parquet",
        },
        {
            "source": "Uni-Mol molecular encoder",
            "molecules_or_samples": "209M conformers in published pretraining; 1,232 local embeddings",
            "targets": "3D denoising and molecular representation pretraining",
            "physics_content": "geometry-aware atom and molecular representations",
            "water_specific": "No",
            "data_class": "pretrained 3D foundation encoder",
            "downloaded": True,
            "usable_for_training": "Yes: frozen embeddings from deterministic RDKit conformers",
            "canonical_url": "https://huggingface.co/dptech/Uni-Mol-Models",
            "doi": "10.1101/2022.11.03.515044",
            "license": "MIT (Uni-Mol repository)",
            "exact_data_fields": "512-dimensional molecular representation",
            "benchmark_overlap": "Pretraining identity unknown; no hydration labels",
            "pretraining_use": "Frozen 3D structure encoder",
            "privileged_supervision": "No",
            "leakage_risk": "No hydration targets; pretraining membership undocumented",
            "local_path": "data/processed/unimol_embeddings.parquet",
        },
        {
            "source": "SPICE",
            "molecules_or_samples": "Large conformer/interaction collection (version-dependent)",
            "targets": "DFT energies, forces, charges, dipoles, bond orders",
            "physics_content": "intramolecular and noncovalent quantum response",
            "water_specific": "No; includes water-containing subsets",
            "data_class": "QM conformers/interactions",
            "downloaded": False,
            "usable_for_training": "Potential foundation pretraining; indirect for hydration",
            "canonical_url": "https://github.com/openmm/spice-dataset",
            "doi": "10.1038/s41597-022-01882-6",
            "license": "CC0-1.0",
            "exact_data_fields": "geometries; energies; forces; multipoles/charges; bond orders",
            "benchmark_overlap": "Unknown; indirect supervision",
            "pretraining_use": "Future 3D physics encoder",
            "privileged_supervision": "Yes",
            "leakage_risk": "Low target leakage",
            "local_path": None,
        },
        {
            "source": "AqSolDB / BigSolDBv2",
            "molecules_or_samples": "9,982 compounds / large multi-solvent solubility corpus",
            "targets": "aqueous or solvent-conditioned logS",
            "physics_content": "solubility (conflates hydration and condensed-phase effects)",
            "water_specific": "AqSolDB yes; BigSolDB no",
            "data_class": "experimental solubility",
            "downloaded": False,
            "usable_for_training": "Only weak auxiliary representation supervision",
            "canonical_url": "https://github.com/mcsorkun/AqSolDB; https://zenodo.org/records/15094979",
            "doi": "10.1038/s41597-019-0151-1",
            "license": "MIT code; dataset/source-specific terms",
            "exact_data_fields": "SMILES; logS; temperature/solvent where available",
            "benchmark_overlap": "Not audited because not selected for direct training",
            "pretraining_use": "Low-priority auxiliary task",
            "privileged_supervision": "No",
            "leakage_risk": "Exact structures must be excluded if used",
            "local_path": None,
        },
        {
            "source": "ChemBERTa-77M-MTR / MoLFormer-XL",
            "molecules_or_samples": "77M / approximately 1.1B pretraining molecules",
            "targets": "masked-token molecular representation",
            "physics_content": "broad chemical structure only",
            "water_specific": "No",
            "data_class": "pretrained foundation encoders",
            "downloaded": True,
            "usable_for_training": "Yes: frozen embeddings then bounded fine-tuning",
            "canonical_url": "https://huggingface.co/DeepChem/ChemBERTa-77M-MTR; https://huggingface.co/ibm-research/MoLFormer-XL-both-10pct",
            "doi": "10.1145/3412841.3440999; 10.1038/s42256-022-00580-7",
            "license": "Model-card-specific; verify before redistribution",
            "exact_data_fields": "pretrained encoder weights/tokenizers",
            "benchmark_overlap": "Pretraining identity unknown; no hydration labels",
            "pretraining_use": "Frozen embeddings",
            "privileged_supervision": "No",
            "leakage_risk": "No target labels, but corpus membership is generally undocumented",
            "local_path": "data/processed/molformer_embeddings.parquet; ChemBERTa checkpoint quarantined by safe-loading policy",
        },
        {
            "source": "OpenADMET Chemprop foundation weights",
            "molecules_or_samples": "1M PubChem pretraining molecules; multiple 34 MB checkpoints",
            "targets": "Jazzy, MiniMol, Mordred, ECFP and other pseudo-descriptors",
            "physics_content": "graph encoders distilled from quantum/property representations",
            "water_specific": "Jazzy checkpoint is hydration-oriented",
            "data_class": "pretrained graph checkpoints",
            "downloaded": True,
            "usable_for_training": "Yes; use as initialization with a leakage caveat for Jazzy",
            "canonical_url": "https://huggingface.co/openadmet/chemprop-foundation-pretraining-weights",
            "doi": "https://arxiv.org/abs/2404.14986",
            "license": "Apache-2.0",
            "exact_data_fields": "Chemprop message-passing state dicts (hidden 2048; depth 6)",
            "benchmark_overlap": "Pretraining corpus identity not published",
            "pretraining_use": "Structure encoder initialization",
            "privileged_supervision": "Jazzy/MiniMol/QM-derived representation targets",
            "leakage_risk": "Jazzy coefficients were fit to Gerber experimental hydration; report separately",
            "local_path": "data/external/huggingface/openadmet__chemprop-foundation-pretraining-weights",
        },
    ]

    expanded = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    abraham = pd.read_parquet(processed / "soluteml_abraham_nonbenchmark.parquet")
    enthalpy = pd.read_parquet(processed / "combisolvh_qm_water_nonbenchmark.parquet")
    acid = pd.read_parquet(processed / "acid_response_nonbenchmark.parquet")
    relative_pairs = pd.read_parquet(processed / "thermodynamic_length_pairs_nonbenchmark.parquet")
    catalog.extend(
        [
            {
                "source": "SoluteML dGsolvDB3 / SoluteDB",
                "molecules_or_samples": (
                    f"{len(expanded):,} strict water hydration values; "
                    f"{len(abraham):,} Abraham-response molecules"
                ),
                "targets": "experimental hydration dG; Abraham E/S/A/B/L",
                "physics_content": (
                    "water partition response and hydrogen-bond/cavity/polarity axes"
                ),
                "water_specific": "Hydration values yes; Abraham axes transferable",
                "data_class": "experimental + empirical physical descriptors",
                "downloaded": True,
                "usable_for_training": "Yes after strict benchmark exclusion",
                "canonical_url": "https://zenodo.org/records/5792296",
                "doi": "10.1021/acs.jcim.1c01103",
                "license": "See Zenodo record and source README",
                "exact_data_fields": (
                    "SMILES; dGsolv; uncertainty/source counts; Abraham E/S/A/B/L"
                ),
                "benchmark_overlap": "0 in processed strict tables",
                "pretraining_use": "Hydration/Abraham multi-task structure encoder",
                "privileged_supervision": "Abraham physical-response heads",
                "leakage_risk": "84 benchmark connectivities removed before aggregation",
                "local_path": "data/processed/expanded_public_hydration_nonbenchmark.parquet",
            },
            {
                "source": "CombiSolvH-QM / SolProp v1.2",
                "molecules_or_samples": f"{len(enthalpy):,} strict water solutes",
                "targets": "COSMO-RS dH; matched dG; T*dS=dH-dG",
                "physics_content": "water-specific thermodynamic response hierarchy",
                "water_specific": "Yes for processed teacher subset",
                "data_class": "QM/continuum thermodynamics",
                "downloaded": True,
                "usable_for_training": "Yes after strict benchmark exclusion",
                "canonical_url": "https://zenodo.org/records/5970538",
                "doi": "10.1021/jacs.2c01768",
                "license": "See Zenodo record and source README",
                "exact_data_fields": "solute/solvent SMILES; dH; matched dG; derived T*dS",
                "benchmark_overlap": "0 in processed strict table",
                "pretraining_use": "Water dG/dH/TdS hierarchy encoder",
                "privileged_supervision": "Enthalpic and entropic response heads",
                "leakage_risk": "All benchmark solutes removed before matching",
                "local_path": "data/processed/combisolvh_qm_water_nonbenchmark.parquet",
            },
            {
                "source": "DISSOLVE2-ANIONS",
                "molecules_or_samples": f"{len(acid):,} strict neutral acid structures",
                "targets": "aqueous pKa; gas-phase acidity; COSMO-RS anion solvation",
                "physics_content": "deprotonation and charged water-response hierarchy",
                "water_specific": "Yes",
                "data_class": "experimental + QM/continuum",
                "downloaded": True,
                "usable_for_training": "Safe targets only; upstream neutral dG quarantined",
                "canonical_url": "https://zenodo.org/records/15604045",
                "doi": "10.5281/zenodo.15604045",
                "license": "See Zenodo record",
                "exact_data_fields": "neutral SMILES; pKa; gas acidity; anion COSMO-RS dG",
                "benchmark_overlap": "0 in processed strict table",
                "pretraining_use": "Acid/anion water-response auxiliary heads",
                "privileged_supervision": "Yes",
                "leakage_risk": (
                    "DirectML neutral dG and derived ion dG excluded from headline features"
                ),
                "local_path": "data/processed/acid_response_nonbenchmark.parquet",
            },
            {
                "source": "Relative-solvation multi-code archive",
                "molecules_or_samples": "9 absolute solutes plus relative transformations",
                "targets": "lambda-resolved gradients for archived runs; multi-code dG",
                "physics_content": "explicit alchemical response under AMBER/CHARMM/GROMACS/SOMD",
                "water_specific": "Yes",
                "data_class": "classical explicit-solvent MD",
                "downloaded": True,
                "usable_for_training": "Too few distinct solutes for a general lambda model",
                "canonical_url": "https://github.com/halx/relative-solvation-inputs",
                "doi": "10.1021/acs.jctc.8b00544",
                "license": "No repository license found",
                "exact_data_fields": "simulation inputs; u_kl; gradients; selected TI curves",
                "benchmark_overlap": "Most small absolute solutes overlap the benchmark",
                "pretraining_use": "Audit/reference only at current molecule count",
                "privileged_supervision": "Lambda-response curves where present",
                "leakage_risk": "Must outer-fold mask overlapping benchmark molecules",
                "local_path": "data/external/repos/relative-solvation-inputs",
            },
            {
                "source": "Thermodynamic-length relative hydration pairs",
                "molecules_or_samples": (
                    f"{len(relative_pairs)} strict pairs / "
                    f"{len(set(relative_pairs.connectivity_key_a) | set(relative_pairs.connectivity_key_b))} molecules"
                ),
                "targets": "vacuum/solvent relative dG and uncertainty",
                "physics_content": "phase-resolved relative alchemical response",
                "water_specific": "Yes",
                "data_class": "classical explicit-solvent MD",
                "downloaded": True,
                "usable_for_training": "Yes as a small pairwise/contrastive auxiliary task",
                "canonical_url": "https://github.com/choderalab/thermo-length-learning",
                "doi": "10.48550/arXiv.1906.08599",
                "license": "No repository license found",
                "exact_data_fields": "pair SMILES; vacuum/solvent dG; phase uncertainty",
                "benchmark_overlap": "0; 32 overlapping pairs removed",
                "pretraining_use": "Relative-response and uncertainty supervision",
                "privileged_supervision": "Yes",
                "leakage_risk": "Both pair endpoints screened by connectivity",
                "local_path": "data/processed/thermodynamic_length_pairs_nonbenchmark.parquet",
            },
            {
                "source": "SAMPL4 AMOEBA polarizable hydration protocols",
                "molecules_or_samples": (
                    f"{len(sampl4_physics)} strict SAMPL4 structures; "
                    "46 with seven triplicate AMOEBA/GAFF protocols"
                ),
                "targets": (
                    "AMOEBA Poltype/polarization-group/multiconformer/basis/water/"
                    "OH-scale dG; GAFF dG; protocol spread and deltas"
                ),
                "physics_content": (
                    "polarizable versus fixed-charge explicit-water alchemical response"
                ),
                "water_specific": "Yes",
                "data_class": "classical polarizable/fixed-charge MD",
                "downloaded": True,
                "usable_for_training": "Yes as benchmark-disjoint privileged targets",
                "canonical_url": "https://zenodo.org/records/35586",
                "doi": "10.5281/zenodo.35586",
                "license": sampl4_metadata["amoeba_license"],
                "exact_data_fields": ";".join(sampl4_physics.columns),
                "benchmark_overlap": "1 connectivity (pyrene) removed; 0 in strict table",
                "pretraining_use": "Polarizable-water response heads",
                "privileged_supervision": "Yes",
                "leakage_risk": "All benchmark connectivities removed before teacher fitting",
                "local_path": "data/processed/sampl4_physics_nonbenchmark.parquet",
            },
            {
                "source": "Bannan et al. partition-coefficient alchemical archive",
                "molecules_or_samples": (
                    f"{len(bannan_curves)} strict of "
                    f"{bannan_metadata['source_rows_before_benchmark_exclusion']} solutes; "
                    "19 water alchemical edges per solute"
                ),
                "targets": (
                    "19 adjacent-state MBAR increments and uncertainties; GAFF/DC totals; "
                    "curve PCs; water/cyclohexane/octanol archive"
                ),
                "physics_content": "explicit-solvent alchemical response curves",
                "water_specific": "Processed teacher is water; archive includes two organic solvents",
                "data_class": "classical explicit-solvent MD",
                "downloaded": True,
                "usable_for_training": "Yes as a small lambda-response teacher",
                "canonical_url": "https://zenodo.org/records/4977926",
                "doi": bannan_metadata["doi"],
                "license": bannan_metadata["license"],
                "exact_data_fields": ";".join(bannan_curves.columns),
                "benchmark_overlap": (
                    f"{bannan_metadata['benchmark_connectivity_overlap_removed']} removed; "
                    "0 in strict table"
                ),
                "pretraining_use": "Alchemical-curve and force-field-delta heads",
                "privileged_supervision": "Yes",
                "leakage_risk": "All benchmark connectivities removed before PCA/teacher fitting",
                "local_path": "data/processed/bannan_alchemical_curve_nonbenchmark.parquet",
            },
            {
                "source": "AMOEBA multi-solvent small-molecule solvation",
                "molecules_or_samples": "21 solutes in toluene/chloroform; 6 in acetonitrile/DMSO",
                "targets": "experimental, triplicate AMOEBA, and triplicate GAFF solvation dG",
                "physics_content": "polarizable/fixed-charge response across four nonaqueous solvents",
                "water_specific": "No",
                "data_class": "experimental + classical polarizable/fixed-charge MD",
                "downloaded": True,
                "usable_for_training": "Potential solvent-conditioned auxiliary task; too small alone",
                "canonical_url": "https://zenodo.org/records/59203",
                "doi": "10.5281/zenodo.59203",
                "license": "CC-BY-4.0",
                "exact_data_fields": "experimental dG; three AMOEBA repeats; three GAFF repeats",
                "benchmark_overlap": "Not used in current water teacher",
                "pretraining_use": "Polarizable-minus-fixed-charge solvent response",
                "privileged_supervision": "Yes",
                "leakage_risk": "Must map names and exclude benchmark solutes before use",
                "local_path": "data/external/zenodo/59203/RESULTS",
            },
            {
                "source": "Replica Exchange with Flexible Timing hydration benchmark",
                "molecules_or_samples": (
                    f"{reft.Molecule.nunique()} molecules / {len(reft)} convergence rows"
                ),
                "targets": "hydration dG and uncertainty versus method and ns/window",
                "physics_content": "sampling/convergence response across REFT, FEP, and Transformato",
                "water_specific": "Yes",
                "data_class": "classical explicit-water alchemical MD",
                "downloaded": True,
                "usable_for_training": "Diagnostic convergence supervision; only seven molecules",
                "canonical_url": "https://zenodo.org/records/20699238",
                "doi": "10.5281/zenodo.20699238",
                "license": "CC-BY-4.0",
                "exact_data_fields": ";".join(reft.columns),
                "benchmark_overlap": "Requires identity resolution before supervised use",
                "pretraining_use": "Sampling-error/convergence head only",
                "privileged_supervision": "Yes",
                "leakage_risk": "Names require structure resolution and benchmark exclusion",
                "local_path": "data/external/zenodo/20699238/hydration_free_enenergies.csv",
            },
            {
                "source": "OpenFE/OpenFF 2.3.0 FreeSolv ASFE archive",
                "molecules_or_samples": (
                    "603 solutes x 3 repeats x solvent/vacuum legs; "
                    f"{len(openfe_diagnostics)} strict nonbenchmark solutes"
                ),
                "targets": (
                    "hydration and leg dG; uncertainties; 14-state MBAR overlap matrices; "
                    "forward/reverse convergence; replica-exchange spectra; NAGL charges"
                ),
                "physics_content": (
                    "10 ns explicit-water alchemical response, convergence, and mixing diagnostics"
                ),
                "water_specific": "Yes",
                "data_class": "classical explicit-water alchemical MD",
                "downloaded": True,
                "usable_for_training": "Yes as benchmark-disjoint privileged supervision",
                "canonical_url": "https://zenodo.org/records/21810272",
                "doi": openfe_metadata["doi"],
                "license": openfe_metadata["license"],
                "exact_data_fields": ";".join(openfe_diagnostics.columns),
                "benchmark_overlap": "83 source records removed; 0 in strict table",
                "pretraining_use": "Classical-response and simulation-quality hierarchy heads",
                "privileged_supervision": "Yes",
                "leakage_risk": (
                    "All benchmark connectivities removed before diagnostic teacher fitting"
                ),
                "local_path": "data/processed/openfe_diagnostics_nonbenchmark.parquet",
            },
            {
                "source": "NQELiq-298",
                "molecules_or_samples": (
                    f"92 molecular liquids; {len(nqeliq)} strict nonbenchmark structures"
                ),
                "targets": (
                    "classical/PIMD-H/PIMD-D density, volume, expansivity, compressibility, "
                    "dielectric constant, dHvap; NQE and isotope shifts"
                ),
                "physics_content": "paired classical and path-integral bulk-liquid response",
                "water_specific": "No; chemically diverse neat molecular liquids",
                "data_class": "classical MD + PIMD",
                "downloaded": True,
                "usable_for_training": "Yes as benchmark-disjoint NQE auxiliary supervision",
                "canonical_url": "https://zenodo.org/records/15236881",
                "doi": nqeliq_metadata["doi"],
                "license": nqeliq_metadata["license"],
                "exact_data_fields": ";".join(nqeliq.columns),
                "benchmark_overlap": (
                    f"{nqeliq_metadata['benchmark_connectivity_overlap_removed']} removed; "
                    "0 in strict table"
                ),
                "pretraining_use": "Narrow classical-to-PIMD/NQE hierarchy heads",
                "privileged_supervision": "Yes",
                "leakage_risk": "All benchmark connectivities removed before teacher fitting",
                "local_path": "data/processed/nqeliq_nonbenchmark.parquet",
            },
            {
                "source": "MolSolv SMD(water)",
                "molecules_or_samples": (
                    f"{molsolv_metadata['records_total']:,} source conformers; "
                    f"{len(molsolv_smd):,} strict sampled structures"
                ),
                "targets": "M06-2X/6-31G* SMD(water) solvation free energy",
                "physics_content": (
                    "large water-specific quantum-continuum solvent-response teacher"
                ),
                "water_specific": "Yes",
                "data_class": "DFT/SMD continuum water",
                "downloaded": True,
                "usable_for_training": "Yes after global benchmark exclusion",
                "canonical_url": "https://zenodo.org/records/7262826",
                "doi": molsolv_metadata["doi"],
                "license": molsolv_metadata["license"],
                "exact_data_fields": ";".join(molsolv_smd.columns),
                "benchmark_overlap": (
                    f"{molsolv_metadata['benchmark_exact_structure_matches_removed']} "
                    "structures removed; 0 in strict table"
                ),
                "pretraining_use": "Water-response encoder pretraining",
                "privileged_supervision": "Yes: approximate water free energy",
                "leakage_risk": (
                    "All benchmark identities removed before sampling, deduplication, and fitting"
                ),
                "local_path": "data/processed/molsolv_smd_water_nonbenchmark.parquet",
            },
            {
                "source": "FreeSolv explicit-solvent phase-space averages",
                "molecules_or_samples": (
                    f"642 source molecules; {len(phase_space)} strict nonbenchmark connectivities"
                ),
                "targets": "Boltzmann-averaged intramolecular distance matrices",
                "physics_content": (
                    "explicit-water trajectory-averaged geometry, conformational response, "
                    "distance/Coulomb spectra"
                ),
                "water_specific": "Yes",
                "data_class": "explicit-water MD phase-space averages",
                "downloaded": True,
                "usable_for_training": "Yes as benchmark-disjoint privileged supervision",
                "canonical_url": "https://zenodo.org/records/6401711",
                "doi": "10.5281/zenodo.6401711",
                "license": "CC-BY-4.0",
                "exact_data_fields": ";".join(phase_space.columns),
                "benchmark_overlap": "80 source records removed; 0 in strict table",
                "pretraining_use": "Structure-to-solvent-averaged-geometry response heads",
                "privileged_supervision": "Yes",
                "leakage_risk": (
                    "Benchmark connectivities removed before response aggregation and fitting"
                ),
                "local_path": "data/processed/phase_space_dynamic_nonbenchmark.parquet",
            },
            {
                "source": "G-NequIP paired gas/SMD solvation benchmark",
                "molecules_or_samples": (
                    "428 energy records / 424 paired geometries; "
                    f"{len(gnequip)} strict nonbenchmark connectivities"
                ),
                "targets": (
                    "QM SMD dG; NNP solvent/geometry/total response; "
                    "solvent-induced geometry changes"
                ),
                "physics_content": "paired gas-to-water implicit-solvent energy hierarchy",
                "water_specific": "Yes",
                "data_class": "DFT/SMD + equivariant NNP",
                "downloaded": True,
                "usable_for_training": "Yes as benchmark-disjoint privileged supervision",
                "canonical_url": "https://zenodo.org/records/20690503",
                "doi": "10.5281/zenodo.20690503",
                "license": "Zenodo record terms; verify before redistribution",
                "exact_data_fields": ";".join(gnequip.columns),
                "benchmark_overlap": "71 source connectivities removed; 0 in strict table",
                "pretraining_use": "Narrow gas→SMD→QM/NNP response hierarchy",
                "privileged_supervision": "Yes",
                "leakage_risk": "All benchmark connectivities removed before teacher fitting",
                "local_path": "data/processed/gnequip_solvation_nonbenchmark.parquet",
            },
            {
                "source": "MLFF hydration free-energy comparison",
                "molecules_or_samples": f"59 source solutes; {len(mlff_hfe)} strict rows",
                "targets": "Organic-MPNICE, GAFF, OPLS4, E-sol, and DFT/PBF hydration dG",
                "physics_content": "matched multi-fidelity force-field hydration hierarchy",
                "water_specific": "Yes",
                "data_class": "MLFF + classical force fields + continuum DFT",
                "downloaded": True,
                "usable_for_training": "Yes, but only 45 strict structures",
                "canonical_url": "https://figshare.com/articles/dataset/31809957",
                "doi": "10.1021/acs.jctc.5c02019",
                "license": "Figshare record terms",
                "exact_data_fields": ";".join(mlff_hfe.columns),
                "benchmark_overlap": "14 source connectivities removed; 0 in strict table",
                "pretraining_use": "Matched MLFF/classical/DFT hierarchy heads",
                "privileged_supervision": "Yes",
                "leakage_risk": "All benchmark connectivities removed before fitting",
                "local_path": "data/processed/mlff_hfe_nonbenchmark.parquet",
            },
            {
                "source": "FreeSolv archived per-window GROMACS energies",
                "molecules_or_samples": "Documented for the 2017 642-molecule campaign",
                "targets": "20-state fep/vdW dH/dlambda and XVG energy time series",
                "physics_content": "full classical alchemical response curves",
                "water_specific": "Yes",
                "data_class": "classical explicit-water alchemical MD",
                "downloaded": False,
                "usable_for_training": "Not currently: documented directory absent from archives",
                "canonical_url": "https://github.com/MobleyLab/FreeSolv/issues/52",
                "doi": "10.5281/zenodo.1161245",
                "license": "FreeSolv CC-BY-4.0 if recovered",
                "exact_data_fields": "XVG dH/dlambda at 20 lambda states",
                "benchmark_overlap": "Would require the same 80-connectivity exclusion",
                "pretraining_use": "Highest-priority lambda-response teacher if recovered",
                "privileged_supervision": "Yes",
                "leakage_risk": "Must exclude benchmark molecules before curve model fitting",
                "local_path": "Unavailable: gromacs_energies omitted from GitHub/Zenodo ZIP",
            },
        ]
    )

    gnnis_metadata = json.loads(
        (processed / "gnnis_static_response_features.metadata.json").read_text()
    )
    lambda_potential_metadata = json.loads(
        (processed / "lambda_potential_static_features.metadata.json").read_text()
    )
    catalog.extend(
        [
            {
                "source": "GNNImplicitSolvent explicit-water mean-force model",
                "molecules_or_samples": (
                    "369,486 molecules / approximately 3.2M conformers; "
                    f"{gnnis_metadata['successful_molecules']:,} local static evaluations"
                ),
                "targets": "explicit-water mean solvation forces and learned solvent energy",
                "physics_content": "force-matched solvent response from explicit-water trajectories",
                "water_specific": "Yes",
                "data_class": "explicit-water MD force teacher + public GNN checkpoint",
                "downloaded": True,
                "usable_for_training": "Yes; checkpoint/static response is benchmark-label independent",
                "canonical_url": "https://github.com/rinikerlab/GNNImplicitSolvent",
                "doi": "10.1039/D4SC02432J; 10.3929/ethz-b-000667722",
                "license": "CC-BY-SA-4.0 data; repository terms apply to code",
                "exact_data_fields": "static energy, force, response, and hidden-state summaries",
                "benchmark_overlap": "No hydration labels used; 85 benchmark structures evaluated statically",
                "pretraining_use": "Explicit-water solvent-force representation",
                "privileged_supervision": "Yes: force/energy response teacher",
                "leakage_risk": "No experimental hydration targets in checkpoint training",
                "local_path": "data/processed/gnnis_static_response_features.parquet",
            },
            {
                "source": "Lambda-aware implicit-solvent force checkpoint",
                "molecules_or_samples": (
                    "280K training configurations reported; "
                    f"{lambda_potential_metadata['successful_molecules']:,} local static evaluations"
                ),
                "targets": "solvent forces plus steric/electrostatic lambda derivatives",
                "physics_content": "learned lambda-resolved solvent response and endpoint energy",
                "water_specific": "Yes",
                "data_class": "explicit-solvent force/derivative distillation checkpoint",
                "downloaded": True,
                "usable_for_training": "Yes; deterministic structure/conformer evaluation",
                "canonical_url": "https://github.com/Popov-Lab-UNC/mlimplicitsolvent-RAW",
                "doi": None,
                "license": "Repository license/provenance; verify before redistribution",
                "exact_data_fields": "steric/electrostatic lambda energies, finite-difference responses, forces",
                "benchmark_overlap": "No experimental hydration target used in the extracted response",
                "pretraining_use": "Lambda-response auxiliary head or static teacher",
                "privileged_supervision": "Yes",
                "leakage_risk": "Raw 280K training identities are not published with the checkpoint",
                "local_path": "data/processed/lambda_potential_static_features.parquet",
            },
            {
                "source": "ConfSolv COSMO-RS conformer solvation",
                "molecules_or_samples": (
                    f"{confsolv_metadata['raw_water_rows']:,} H2O conformer rows; "
                    f"{len(confsolv_water):,} strict neutral connectivities"
                ),
                "targets": "Boltzmann-ensemble hydration plus conformer-resolved COSMO-RS response",
                "physics_content": "gas/solution conformer corrections and water response distributions",
                "water_specific": "Yes (H2O slice extracted from 41-solvent archive)",
                "data_class": "DFT/COSMO-RS conformer ensemble",
                "downloaded": True,
                "usable_for_training": "Yes; structure mapping and benchmark exclusion complete",
                "canonical_url": "https://zenodo.org/records/8292520",
                "doi": "10.5281/zenodo.8292520",
                "license": "CC-BY-4.0",
                "exact_data_fields": ";".join(confsolv_water.columns),
                "benchmark_overlap": (
                    f"{confsolv_metadata['benchmark_connectivity_overlaps_removed']} "
                    "connectivities removed globally"
                ),
                "pretraining_use": "Conformer-aware water-response hierarchy",
                "privileged_supervision": "Yes: structure-only surrogate targets",
                "leakage_risk": "Resolved by connectivity exclusion; unverified ID variants omitted",
                "local_path": "data/processed/confsolv_water_nonbenchmark.parquet",
            },
            {
                "source": "ReSolv implicit-solvent free-energy potential",
                "molecules_or_samples": "583 paired vacuum/water trajectories; 389 FreeSolv training molecules",
                "targets": "experimental hydration free energy through differentiable BAR training",
                "physics_content": "vacuum/water conformer distributions and learned solvent potential",
                "water_specific": "Yes",
                "data_class": "ML potential + trajectory reweighting",
                "downloaded": True,
                "usable_for_training": "No for headline benchmark; checkpoint is target-contaminated",
                "canonical_url": "https://github.com/tummfm/ReSolv",
                "doi": "10.1063/5.0235189",
                "license": "Repository license; Git LFS budget currently exhausted",
                "exact_data_fields": "paired precomputed trajectories and trained U_vac/U_wat checkpoints",
                "benchmark_overlap": "FreeSolv ancestry includes many ARROW benchmark molecules",
                "pretraining_use": "Diagnostic only unless retrained after solute exclusion",
                "privileged_supervision": "Potentially, after leakage-safe retraining",
                "leakage_risk": "Published U_wat checkpoint was trained on experimental FreeSolv labels",
                "local_path": "data/external/repos/ReSolv",
            },
            {
                "source": "Freecurve all-branch Git physics audit",
                "molecules_or_samples": "89 branches / 7,513 physics-related path matches",
                "targets": "TI/BAR/dHdl/PIMD/energy output discovery",
                "physics_content": "historical Git-tree audit across priority repositories",
                "water_specific": "Mixed",
                "data_class": "repository provenance audit",
                "downloaded": True,
                "usable_for_training": "Only identified molecule-resolved outputs",
                "canonical_url": "private Freecurve organization",
                "doi": None,
                "license": "Internal Freecurve",
                "exact_data_fields": "repository, branch, commit, path, blob size",
                "benchmark_overlap": "Toluene plus two isobutylbenzene endpoints; one unidentified curve",
                "pretraining_use": "Data-availability proof, not a model target by itself",
                "privileged_supervision": "Sparse identified outputs only",
                "leakage_risk": "Unidentified campaign quarantined",
                "local_path": "data/catalog/freecurve_all_branch_physics_paths.tsv",
            },
        ]
    )

    catalog_frame = pd.DataFrame(catalog)
    # This field intentionally mixes exact integer counts with audit-status
    # text for sources that have not been materialized molecule-by-molecule.
    catalog_frame["benchmark_overlap"] = catalog_frame["benchmark_overlap"].astype("string")
    catalog_frame["retrieved_utc"] = pd.Timestamp.utcnow().isoformat()
    catalog_frame.to_parquet(catalog_dir / "global_data_catalog.parquet", index=False)

    summary = {
        "combisolv_qm_rows": len(combi_qm),
        "combisolv_qm_unique_solutes": int(combi_qm.solute_connectivity_key.nunique()),
        "combisolv_qm_unique_solvents": int(combi_qm.solvent_connectivity_key.nunique()),
        "combisolv_qm_water_rows": int(combi_qm.solvent_canonical_smiles.eq("O").sum()),
        "combisolv_qm_benchmark_overlap_rows": int(combi_qm.benchmark_solute_overlap.sum()),
        "combisolv_qm_strict_rows": len(strict_qm),
        "combisolv_qm_strict_water_rows": len(water_qm),
        "combisolv_exp_rows": len(combi_exp),
        "combisolv_exp_benchmark_overlap_rows": int(combi_exp.benchmark_solute_overlap.sum()),
        "g4mp2_rows": len(g4),
        "g4mp2_benchmark_overlap_rows": int(g4.benchmark_solute_overlap.sum()),
        "g4mp2_strict_rows": len(strict_g4),
        "guthrie_curated_rows": len(guthrie),
        "guthrie_benchmark_overlap_rows": int(guthrie.benchmark_solute_overlap.sum()),
        "sampl1_rows": len(sampl),
        "sampl1_benchmark_overlap_rows": int(sampl.benchmark_solute_overlap.sum()),
        "solvatum_solutes": solvatum_molecules,
        "solvatum_pairs": solvatum_pairs,
        "thermodynamic_length_molecules": thermo_molecules,
        "thermodynamic_length_result_directories": thermo_pairs,
        "probe_completed": probe_count,
        "expanded_hydration_strict_rows": len(expanded),
        "abraham_strict_rows": len(abraham),
        "combisolvh_qm_water_strict_rows": len(enthalpy),
        "acid_response_strict_rows": len(acid),
        "relative_alchemical_safe_pairs": len(relative_pairs),
        "des370k_water_geometries": des_metadata["water_dimer_geometries"],
        "des370k_water_strict_solutes": len(des_water),
        "des370k_benchmark_overlap_removed": des_metadata["benchmark_connectivity_overlap_removed"],
        "sampl4_physics_strict_rows": len(sampl4_physics),
        "sampl4_amoeba_label_rows": int(sampl4_physics.amoeba_consensus.notna().sum()),
        "bannan_alchemical_strict_rows": len(bannan_curves),
        "bannan_alchemical_edges_per_molecule": bannan_metadata["raw_curve_edges_per_molecule"],
        "reft_convergence_rows": len(reft),
        "reft_convergence_molecules": int(reft.Molecule.nunique()),
        "openfe_diagnostics_strict_rows": len(openfe_diagnostics),
        "openfe_diagnostics_unit_results": openfe_metadata["protocol_unit_results"],
        "nqeliq_strict_rows": len(nqeliq),
        "nqeliq_benchmark_overlap_removed": nqeliq_metadata[
            "benchmark_connectivity_overlap_removed"
        ],
        "molsolv_smd_source_records": molsolv_metadata["records_total"],
        "molsolv_smd_strict_structures": len(molsolv_smd),
        "molsolv_smd_benchmark_removed": molsolv_metadata[
            "benchmark_exact_structure_matches_removed"
        ],
        "phase_space_strict_rows": len(phase_space),
        "phase_space_targets": phase_space_metadata["target_count"],
        "gnequip_strict_rows": len(gnequip),
        "gnequip_complete_rows": gnequip_metadata["complete_source_connectivities"],
        "mlff_hfe_strict_rows": len(mlff_hfe),
        "mlff_hfe_benchmark_removed": mlff_hfe_metadata["benchmark_connectivities_removed"],
        "confsolv_water_raw_rows": confsolv_metadata["raw_water_rows"],
        "confsolv_water_strict_connectivities": len(confsolv_water),
        "confsolv_benchmark_removed": confsolv_metadata["benchmark_connectivity_overlaps_removed"],
        "source_commits": {
            path.name: git_commit(path)
            for path in (ROOT / "data/external/repos").iterdir()
            if path.is_dir()
        },
    }
    (catalog_dir / "global_data_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
