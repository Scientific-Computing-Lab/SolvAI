#!/usr/bin/env python3
"""Install the frozen Penpot Figure 1 exports into paper-facing locations."""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PENPOT = ROOT / "figures" / "penpot"
EXPORTS = PENPOT / "exports"
MAIN = ROOT / "figures" / "main"
PAPER = ROOT / "paper" / "figures" / "main"
ALTERNATIVES = ROOT / "figures" / "alternatives"


def verify_frozen_exports() -> None:
    """Require portable, already-normalized Penpot exports.

    Normalization is an explicit authoring step because PDF/PNG renderers can
    produce byte-level differences across platforms.  Paper reproduction copies
    the frozen exports verbatim instead of rerendering them in CI.
    """
    for svg_path in sorted(EXPORTS.glob("fig1_variant_*.svg")):
        svg = svg_path.read_text(encoding="utf-8")
        if 'width="180mm"' not in svg or 'height="112mm"' not in svg:
            raise RuntimeError(
                f"Penpot export is not normalized to publication size: {svg_path}"
            )
        if "http://localhost:9001" in svg:
            raise RuntimeError(f"Unresolved local Penpot URL in {svg_path}")


def copy_triplet(stem: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf", "png"):
        source = EXPORTS / f"{stem}.{suffix}"
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, destination / source.name)


def main() -> None:
    verify_frozen_exports()

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
