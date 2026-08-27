# Editorial summary

Accurate hydration free energies usually require explicit-solvent sampling, and
nuclear quantum effects make the highest-fidelity calculations especially
costly. SolvAI changes when that cost is paid. It learns compact,
structure-predictable solvent-response coordinates from large
benchmark-disjoint quantum-continuum, conformational and alchemical datasets,
then uses those predicted responses to estimate hydration free energy directly
from a SMILES string. On the same chemically diverse 85-solute reference set
used to establish the ARROW/PIMD8 result, SolvAI attains a strict five-fold OOF
MAE of 0.197 kcal/mol, compared with 0.205 kcal/mol for PIMD8 and 0.239 kcal/mol
for the previous structure-only baseline. Five independent partitions average
0.204 ± 0.005 kcal/mol, so the work claims PIMD8-level accuracy rather than a
general sub-0.20 threshold. The broader advance is a practical relationship
between simulation and AI: expensive physical calculations become reusable
training supervision rather than mandatory per-molecule computations.
