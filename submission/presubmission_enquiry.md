# Presubmission enquiry

**Proposed title:** Structure-predicted solvent responses enable simulation-free
hydration free-energy prediction

We ask whether costly physical calculations can define response coordinates that are
learned once and reused across molecules. SolvAI predicts 15 solvent-response priors
from structure and combines them with molecular descriptors in an experimentally
supervised hydration endpoint. PIMD is not a retained teacher: it is the high-fidelity
reference against which structure-only deployment is compared. In a preregistered,
fully matched five-fold analysis on the 85-solute ARROW set, aligned response priors
reduce MAE from 0.303 to 0.202 kcal mol⁻¹ (95% paired-bootstrap interval for the
change, −0.215 to −0.020), comparable to ARROW/PIMD8 at 0.205 kcal mol⁻¹. Shuffled
priors abolish the gain. The advantage persists when related endpoint-labelled
chemistry is removed globally and when all ARROW experimental labels are withheld
from training. Conversely, sparse high-fidelity response targets fail because they
cannot yet be inferred accurately enough from structure. The work therefore provides
both a positive demonstration and a mechanistic boundary for using physical
calculation as reusable supervision rather than per-query computation.
