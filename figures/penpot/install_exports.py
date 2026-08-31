#!/usr/bin/env python3
"""Install the frozen Penpot Figure 1 exports into paper-facing locations."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PENPOT = ROOT / "figures" / "penpot"
EXPORTS = PENPOT / "exports"
MAIN = ROOT / "figures" / "main"
PAPER = ROOT / "paper" / "figures" / "main"
ALTERNATIVES = ROOT / "figures" / "alternatives"


def copy_triplet(stem: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf", "png"):
        source = EXPORTS / f"{stem}.{suffix}"
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, destination / source.name)


def main() -> None:
    subprocess.run([sys.executable, str(PENPOT / "normalize_exports.py")], check=True)

    selected = "fig1_variant_C_balanced"
    copy_triplet(selected, MAIN)
    for suffix in ("svg", "pdf", "png"):
        source = EXPORTS / f"{selected}.{suffix}"
        shutil.copy2(source, MAIN / f"fig1_concept.{suffix}")
        PAPER.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, PAPER / f"fig1_concept.{suffix}")

    for alternative in ("fig1_variant_A_minimal", "fig1_variant_B_molecular"):
        copy_triplet(alternative, ALTERNATIVES)

    print("Installed the frozen Penpot Figure 1 exports (variant C selected).")


if __name__ == "__main__":
    main()
