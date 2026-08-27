"""Console entry point."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from .inference import predict_smiles


def parse_smiles(argv: Sequence[str] | None = None) -> list[str]:
    """Parse the public CLI, accepting both ``solvai predict X`` and ``solvai X``."""
    parser = argparse.ArgumentParser(
        description="Predict hydration free energy from SMILES without simulation."
    )
    parser.add_argument("arguments", nargs="+")
    arguments = parser.parse_args(argv).arguments
    if arguments[0] == "predict":
        arguments = arguments[1:]
    if not arguments:
        parser.error("predict requires at least one SMILES string")
    return arguments


def main() -> None:
    smiles_values = parse_smiles()
    predictions, spread = predict_smiles(smiles_values)
    for smiles, value, uncertainty in zip(smiles_values, predictions, spread, strict=True):
        print(f"{smiles}\t{value:.6f}\t{uncertainty:.6f}")
