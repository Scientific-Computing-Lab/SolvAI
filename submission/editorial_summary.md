# Editorial summary

High-accuracy hydration calculations resolve how a solute reorganizes water, often
through costly sampling repeated for every molecule. SolvAI instead learns 15 compact
solvent-response coordinates from external physical calculations and predicts them
from molecular structure. A separate experimentally supervised model combines these
predicted responses with ordinary molecular descriptors; deployment begins from a
SMILES string and invokes no simulation. In a preregistered matched comparison on the
85-solute ARROW reference set, the response priors lower out-of-fold mean absolute
error from 0.303 to 0.202 kcal mol⁻¹, comparable to the reconstructed 0.205 kcal mol⁻¹
ARROW/PIMD8 result. Shuffled priors do not help, while the advantage persists under
global chemical-separation controls and without ARROW labels in training. PIMD is the
high-fidelity comparator, not a retained teacher. The broader result is a controlled
demonstration that aligned, structure-learnable physical responses can be reused as
molecular supervision instead of recomputed for each query.
