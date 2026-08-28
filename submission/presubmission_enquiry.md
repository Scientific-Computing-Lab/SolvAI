# Presubmission enquiry

**Proposed title:** Structure-predicted solvent responses enable simulation-free
hydration free-energy prediction

Prior studies have used physical descriptors at prediction time and have learned
quantum descriptors from structure. We ask the next solvation-specific question:
whether complementary responses learned from calculated, empirical and corrected
solvation sources can become reusable molecular supervision. SolvAI predicts 15 such
priors from structure and combines them with molecular descriptors in an
experimentally supervised hydration endpoint. PIMD is not a retained teacher: it is
the high-fidelity reference against which structure-only deployment is compared. In a
preregistered, fully matched five-fold analysis on the 85-solute ARROW set, aligned
response priors reduce MAE from 0.303 to 0.202 kcal mol⁻¹ (95% paired-bootstrap
interval for the change, −0.215 to −0.020), comparable to ARROW/PIMD8 at 0.205 kcal
mol⁻¹. Shuffled priors abolish the gain. The advantage persists when related
endpoint-labelled chemistry is removed globally and when all ARROW experimental
labels are withheld from training. Conversely, sparse high-fidelity response targets
fail because they cannot yet be inferred accurately enough from structure. The work
therefore provides both a positive demonstration and a mechanistic boundary for
converting physical-response information into reusable supervision rather than
per-query computation.
