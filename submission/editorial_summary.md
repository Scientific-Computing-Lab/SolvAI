# Editorial summary

High-accuracy hydration calculations resolve how a solute reorganizes water, often
through costly sampling repeated for every molecule. SolvAI consolidates calculated,
empirical and corrected solvent responses into 15 compact coordinates and learns to
predict them from molecular structure. A separate experimentally supervised model
combines these priors with ordinary molecular descriptors; deployment begins from a
SMILES string and invokes no simulation. In a preregistered matched comparison on the
85-solute ARROW reference set, the response layer lowers out-of-fold mean absolute
error from 0.303 to 0.202 kcal mol⁻¹, comparable to the reconstructed 0.205 kcal mol⁻¹
ARROW/PIMD8 result. Shuffled priors do not help, while the advantage persists under
global chemical-separation controls and without ARROW labels in training. PIMD is the
high-fidelity comparator, not a retained teacher. In a prospectively frozen external
molecule-disjoint cohort, the matched response layer lowers MAE from 1.532 to 1.153
kcal mol⁻¹; the advantage remains in 97 molecules absent from every response source,
although absolute accuracy is lower. The broader result is a controlled demonstration
that structure-learnable solvent responses can convert costly physical information
into reusable molecular supervision without recomputation for each query.
