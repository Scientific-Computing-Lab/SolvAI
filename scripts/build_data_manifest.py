#!/usr/bin/env python3
"""Freeze hashes and redistribution decisions for final-model training sources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


SOURCES = [
    {
        "source": "MolSolv SMD(water)",
        "scientific_role": "water-specific quantum-continuum response teacher",
        "processed_file": "molsolv_smd_water_nonbenchmark.parquet",
        "original_records": 1_729_545,
        "filtered_records": 350_391,
        "benchmark_overlaps_removed": 82,
        "version": "Zenodo record 7262826",
        "url": "https://zenodo.org/records/7262826",
        "doi": "10.5281/zenodo.7262826",
        "license": "CC-BY-4.0",
        "redistributed": True,
        "reproduction_command": "python scripts/download_data.py molsolv",
    },
    {
        "source": "ConfSolv H2O",
        "scientific_role": "water conformer and solvent-response hierarchy",
        "processed_file": "confsolv_water_nonbenchmark.parquet",
        "original_records": 5_392_567,
        "filtered_records": 39_878,
        "benchmark_overlaps_removed": 13,
        "version": "Zenodo record 8292520",
        "url": "https://zenodo.org/records/8292520",
        "doi": "10.5281/zenodo.8292520",
        "license": "CC-BY-4.0",
        "redistributed": True,
        "reproduction_command": "python scripts/download_data.py confsolv",
    },
    {
        "source": "SoluteML Abraham",
        "scientific_role": "empirical solute-response axes",
        "processed_file": "soluteml_abraham_nonbenchmark.parquet",
        "original_records": None,
        "filtered_records": 8_098,
        "benchmark_overlaps_removed": 84,
        "version": "Zenodo record 5792296",
        "url": "https://zenodo.org/records/5792296",
        "doi": "10.1021/acs.jcim.1c01103",
        "license": "CC-BY-4.0",
        "redistributed": True,
        "reproduction_command": "python scripts/download_data.py soluteml",
    },
    {
        "source": "OpenFE/OpenFF 2.3.0 FreeSolv ASFE",
        "scientific_role": "explicit-water alchemical response teacher",
        "processed_file": "openff_alchemical_nonbenchmark.parquet",
        "original_records": 603,
        "filtered_records": 520,
        "benchmark_overlaps_removed": 83,
        "version": "Zenodo record 21810272",
        "url": "https://zenodo.org/records/21810272",
        "doi": "10.5281/zenodo.21810272",
        "license": "CC-BY-4.0",
        "redistributed": True,
        "reproduction_command": "python scripts/download_data.py openff",
    },
    {
        "source": "GBn2 / GNNImplicitSolvent",
        "scientific_role": "implicit and learned explicit-water response teacher",
        "processed_file": "implicit_solvent_nonbenchmark.parquet",
        "original_records": 550,
        "filtered_records": 550,
        "benchmark_overlaps_removed": 0,
        "version": "repository state recorded in global catalog",
        "url": "https://github.com/rinikerlab/GNNImplicitSolvent",
        "doi": "10.1039/D4SC02432J",
        "license": "CC-BY-SA-4.0 data",
        "redistributed": True,
        "reproduction_command": "python scripts/download_data.py implicit",
    },
    {
        "source": "CombiSolv-QM water",
        "scientific_role": "COSMOtherm water-response teacher",
        "processed_file": "combisolv_qm_water_nonbenchmark.parquet",
        "original_records": 3_988,
        "filtered_records": 3_963,
        "benchmark_overlaps_removed": 25,
        "version": "Chemical Engineering Journal supplementary file",
        "url": "https://ars.els-cdn.com/content/image/1-s2.0-S1385894721008925-mmc1.txt",
        "doi": "10.1016/j.cej.2021.129307",
        "license": "publisher supplementary terms; no standalone data license",
        "redistributed": False,
        "reproduction_command": "python scripts/download_data.py combisolv-qm --accept-source-terms",
        "processed_bytes": 283_353,
        "processed_sha256": "6fbaceff441b62c8399aae2dbd3b6cfbdfa10b86fc47cee8d60001c73da24832",
    },
]


def main() -> None:
    processed = ROOT / "data/processed"
    records = []
    for source in SOURCES:
        record = dict(source)
        path = processed / str(source["processed_file"])
        if path.exists():
            record["processed_bytes"] = path.stat().st_size
            record["processed_sha256"] = sha256(path)
        elif source["redistributed"]:
            raise FileNotFoundError(path)
        elif "processed_sha256" not in record:
            raise AssertionError(f"Missing frozen hash for non-redistributed source: {path}")
        records.append(record)
    output = ROOT / "data/manifests/training_source_manifest.json"
    output.write_text(json.dumps({"schema_version": 1, "sources": records}, indent=2) + "\n")
    pd.DataFrame(records).to_csv(ROOT / "data/manifests/training_source_manifest.csv", index=False)
    print(f"Wrote {len(records)} source records to {output}")


if __name__ == "__main__":
    main()
