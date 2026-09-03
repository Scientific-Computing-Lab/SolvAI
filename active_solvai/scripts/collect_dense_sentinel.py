"""Collect dense sentinel energy outputs into a versioned response table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from active_solvai.ledger import sha256
from active_solvai.probes import contiguous_block_sem, read_energy

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"


def find_energy(row: dict[str, object]) -> Path:
    if bool(row["existing_observation"]):
        return Path(str(row["energy_file"]))
    candidates = sorted((Path(str(row["case_directory"])) / "output").glob("*__SYSTEM.ene"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected one SYSTEM energy for {row['molecule_name']} lambda={row['lambda']}; "
            f"found {len(candidates)}"
        )
    return candidates[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ACTIVE_ROOT / "simulations/dense_pimd2")
    parser.add_argument("--role", choices=("calibration", "prospective", "all"), default="all")
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    manifest = pd.read_csv(args.root / "manifest.csv")
    if args.role != "all":
        manifest = manifest.loc[manifest.role.eq(args.role)]
    records: list[dict[str, object]] = []
    missing: list[str] = []
    for row in manifest.to_dict("records"):
        try:
            energy = find_energy(row)
        except FileNotFoundError as error:
            missing.append(str(error))
            continue
        frame = read_energy(energy)
        values = frame["dHdL"].dropna().to_numpy(dtype=float)
        sem, blocks = contiguous_block_sem(values)
        records.append(
            {
                "campaign": row["campaign"],
                "role": row["role"],
                "molecule_name": row["molecule_name"],
                "molecule_id": row["molecule_id"],
                "canonical_smiles": row["canonical_smiles"],
                "functional_group_family": row["functional_group_family"],
                "lambda": float(row["lambda"]),
                "ti_point": int(row["ti_point"]),
                "existing_observation": bool(row["existing_observation"]),
                "frames": len(values),
                "mean_dhdl_kcal_mol": float(values.mean()),
                "sd_dhdl_kcal_mol": float(values.std(ddof=1)),
                "five_block_sem_kcal_mol": sem,
                "block_count": blocks,
                "mean_temperature_kelvin": float(frame.Temp.mean()),
                "mean_density_g_cm3": float(frame.Density.mean()),
                "energy_file": str(energy.resolve()),
                "energy_sha256": sha256(energy),
            }
        )
    if missing and not args.allow_incomplete:
        raise RuntimeError("\n".join(missing))
    output_dir = ACTIVE_ROOT / "results/phase2"
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = args.role
    frame = pd.DataFrame(records).sort_values(["role", "molecule_name", "lambda"])
    output = output_dir / f"dense_responses_{suffix}.parquet"
    frame.to_parquet(output, index=False)
    summary = {
        "role": args.role,
        "rows": len(frame),
        "molecules": int(frame.molecule_id.nunique()) if len(frame) else 0,
        "complete_molecules": int((frame.groupby("molecule_id").size() == 15).sum())
        if len(frame)
        else 0,
        "missing": missing,
        "output": str(output),
        "output_sha256": sha256(output),
    }
    (output_dir / f"dense_responses_{suffix}_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
