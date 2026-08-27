#!/usr/bin/env python3
"""Download versioned training sources without placing them under Git control."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "data/external"

ZENODO = {
    "molsolv": ("7262826", None),
    "confsolv": ("8292520", None),
    "soluteml": ("5792296", "Solvation_data-1.0.0.zip"),
    "openff": ("21810272", None),
}
COMBISOLV_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S1385894721008925-mmc1.txt"
COMBISOLV_EXP_URL = (
    "https://raw.githubusercontent.com/su-group/SolvBERT/"
    "77706db56a15403a93a8b9b3083bd385976e1ba0/"
    "solv-bert/data/CombiSolv-Exp-8780.csv"
)
FREESOLV_URL = (
    "https://raw.githubusercontent.com/MobleyLab/FreeSolv/"
    "6c7d19b4b565537365ffd22006aa2cd4643200c6/database.json"
)


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def retrieve(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        print(f"Present: {destination}")
        return
    print(f"Downloading {url} -> {destination}")
    urllib.request.urlretrieve(url, destination)


def download_zenodo(name: str) -> None:
    record, selected = ZENODO[name]
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}") as response:
        metadata = json.load(response)
    destination = EXTERNAL / "zenodo" / record
    for item in metadata["files"]:
        if selected and item["key"] != selected:
            continue
        path = destination / item["key"]
        retrieve(item["links"]["self"], path)
        algorithm, expected = item["checksum"].split(":", 1)
        observed = digest(path, algorithm)
        if observed != expected:
            raise AssertionError(f"Checksum mismatch for {path}")
        print(f"Verified {algorithm}: {observed}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source",
        choices=[*ZENODO, "combisolv-qm", "combisolv-exp", "freesolv", "implicit"],
    )
    parser.add_argument(
        "--accept-source-terms",
        action="store_true",
        help="Required for publisher supplements without a standalone data licence.",
    )
    args = parser.parse_args()
    if args.source in ZENODO:
        download_zenodo(args.source)
    elif args.source in {"combisolv-qm", "combisolv-exp"}:
        if not args.accept_source_terms:
            doi = (
                "10.1016/j.cej.2021.129307"
                if args.source == "combisolv-qm"
                else "10.1039/D2DD00107A"
            )
            raise SystemExit(f"Review DOI {doi} and rerun with --accept-source-terms.")
        if args.source == "combisolv-qm":
            retrieve(COMBISOLV_URL, EXTERNAL / "combisolv-qm" / "CombiSolv-QM.txt")
        else:
            retrieve(
                COMBISOLV_EXP_URL,
                ROOT / "data/raw/public/CombiSolv-Exp-8780.csv",
            )
    elif args.source == "freesolv":
        retrieve(FREESOLV_URL, ROOT / "data/raw/public/freesolv_database.json")
    else:
        raise SystemExit(
            "Clone https://github.com/rinikerlab/GNNImplicitSolvent at the "
            "version recorded in data/manifests/global_data_catalog.parquet."
        )


if __name__ == "__main__":
    main()
