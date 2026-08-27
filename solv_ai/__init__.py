"""SolvAI: structure-only hydration free-energy prediction."""

__version__ = "1.0.0"

from .inference import predict_smiles

__all__ = ["predict_smiles"]
