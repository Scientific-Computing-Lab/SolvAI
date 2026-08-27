"""Console entry point."""

from __future__ import annotations

import argparse

from .inference import predict_smiles


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict hydration free energy from SMILES without simulation."
    )
    parser.add_argument("smiles", nargs="+")
    args = parser.parse_args()
    predictions, spread = predict_smiles(args.smiles)
    for smiles, value, uncertainty in zip(args.smiles, predictions, spread, strict=True):
        print(f"{smiles}\t{value:.6f}\t{uncertainty:.6f}")
