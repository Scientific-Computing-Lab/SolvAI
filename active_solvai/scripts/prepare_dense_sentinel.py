"""Prepare immutable same-Hamiltonian dense PIMD2 sentinel configurations."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

RELEASE_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = RELEASE_ROOT.parents[1]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lambda_slug(value: float) -> str:
    return f"lambda_{value:.2f}".replace(".", "p")


def find_system_energy(output: Path) -> Path:
    candidates = sorted(output.glob("*__SYSTEM.ene"))
    if len(candidates) != 1:
        raise ValueError(f"Expected one SYSTEM energy file in {output}, found {len(candidates)}")
    return candidates[0]


def replace_paths(root: ET.Element, old_case: Path, new_case: Path) -> None:
    for node in root.iter():
        if node.text and str(old_case) in node.text:
            node.text = node.text.replace(str(old_case), str(new_case))
        if "OutputFolder" in node.attrib:
            node.attrib["OutputFolder"] = str(new_case / "output")


def remove_trajectory_outputs(root: ET.Element) -> None:
    for file_output in root.findall(".//FileOutput"):
        for output in list(file_output):
            if output.tag == "Output" and output.attrib.get("DataType") == "TRR":
                file_output.remove(output)


def prepare_new_case(source_case: Path, destination: Path, ti_point: int) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "output").mkdir(exist_ok=True)
    shutil.copyfile(source_case / "solute.hin", destination / "solute.hin")
    tree = ET.parse(source_case / "conf.xml")
    root = tree.getroot()
    replace_paths(root, source_case, destination)
    node = root.find(".//Task[@Type='SetTIPoint']/Settings/Param[@Title='TIPoint']")
    if node is None:
        raise ValueError(f"Missing TIPoint in {source_case / 'conf.xml'}")
    node.text = str(ti_point)
    remove_trajectory_outputs(root)
    ET.indent(tree, space="  ")
    config = destination / "conf.xml"
    tree.write(config, encoding="unicode", xml_declaration=True)
    return config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, default=ACTIVE_ROOT / "configs/dense_sentinel_v1.json"
    )
    parser.add_argument("--output-root", type=Path, default=ACTIVE_ROOT / "simulations/dense_pimd2")
    args = parser.parse_args()
    specification = json.loads(args.config.read_text())
    identity = pd.read_csv(ACTIVE_ROOT / "data/identity/probe_identity_manifest.csv")
    identity = identity.drop_duplicates("molecule_name").set_index("molecule_name")
    roles = {name: "calibration" for name in specification["calibration_molecules"]} | {
        name: "prospective" for name in specification["prospective_molecules"]
    }
    lambda_grid = specification["lambda_grid"]
    existing = specification["existing_lambda_values"]
    campaign_for_lambda = {
        0.1: "pimd2_lambda01_5ps",
        0.5: "pimd2_lambda05_5ps",
        0.9: "pimd2_lambda09_5ps",
    }
    records: list[dict[str, object]] = []
    for name, role in roles.items():
        source_case = (
            WORKSPACE_ROOT
            / "results/physics_probes/pimd2_lambda01_5ps"
            / name.lower().replace(" ", "_")
        )
        if not source_case.is_dir():
            raise FileNotFoundError(source_case)
        for ti_point, lambda_value in enumerate(lambda_grid):
            case = (
                args.output_root / role / name.lower().replace(" ", "_") / lambda_slug(lambda_value)
            )
            is_existing = any(abs(lambda_value - value) < 1e-12 for value in existing)
            if is_existing:
                campaign = campaign_for_lambda[float(lambda_value)]
                historical_case = (
                    WORKSPACE_ROOT / "results/physics_probes" / campaign / source_case.name
                )
                energy = find_system_energy(historical_case / "output")
                config_path = historical_case / "conf.xml"
                case_path = historical_case
            else:
                config_path = prepare_new_case(source_case, case, ti_point)
                energy = case / "output" / "PENDING_SYSTEM.ene"
                case_path = case
            row = identity.loc[name]
            records.append(
                {
                    "campaign": specification["campaign"],
                    "role": role,
                    "molecule_name": name,
                    "molecule_id": row.molecule_id,
                    "canonical_smiles": row.canonical_smiles,
                    "functional_group_family": row.functional_group_family,
                    "lambda": lambda_value,
                    "ti_point": ti_point,
                    "existing_observation": is_existing,
                    "case_directory": str(case_path.resolve()),
                    "config": str(config_path.resolve()),
                    "config_sha256": sha256(config_path),
                    "energy_file": str(energy.resolve()),
                    "production_ps": specification["protocol"]["production_ps_per_window"],
                    "beads": specification["protocol"]["beads"],
                }
            )
    manifest = pd.DataFrame(records).sort_values(["role", "molecule_name", "lambda"])
    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output_root / "manifest.csv", index=False)
    print(
        f"prepared {len(manifest)} total rows; "
        f"{int((~manifest.existing_observation).sum())} new windows"
    )


if __name__ == "__main__":
    main()
