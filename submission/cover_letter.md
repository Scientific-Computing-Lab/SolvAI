# Cover letter

Dear Editors,

We submit “Structure-predicted solvent responses enable simulation-free hydration
free-energy prediction” for consideration as an Article in *Nature Communications*.

High-accuracy solvation calculations repeatedly pay the cost of resolving solvent
response for each molecule. We test a different use of physical computation: learn
compact response coordinates from large external calculations, predict those
coordinates from structure, and expose them to a separately supervised experimental
endpoint. The resulting SolvAI system accepts a SMILES string and runs no molecular
dynamics, path-integral dynamics or probe calculation.

The paper’s central evidence is a preregistered matched comparison. With identical
experimental labels, molecular descriptors, endpoint architecture, weights, folds
and seeds, adding 15 molecule-aligned response priors lowers five-fold out-of-fold MAE
from 0.303 to 0.202 kcal mol⁻¹ on the 85-solute ARROW reference set (paired change
−0.101; 95% bootstrap interval −0.215 to −0.020). The response advantage disappears
when the priors are shuffled, persists across five complete partitions, survives
family, scaffold, molecular-cluster and nearest-neighbour exclusions applied to the
entire supervised pool, and remains when no ARROW experimental labels are used for
training. The final point estimate is comparable to the reconstructed ARROW/PIMD8
error of 0.205 kcal mol⁻¹. PIMD is an accuracy comparator, not a retained teacher.

The work also identifies the boundary of the idea. Compact, learnable solvent-response
coordinates transfer; sparse high-fidelity alchemical responses currently do not,
because their structure-to-response errors remain too large. We believe this
distinction makes the study relevant beyond hydration prediction: it establishes a
controlled route by which physical calculations can become reusable molecular
supervision rather than per-query computation.

The manuscript is accompanied by molecule-level predictions, the preregistration,
identity and similarity audits, all negative controls, frozen artifacts, source code
and a complete reproducibility package.

Sincerely,

Gal Oren and Michael Levitt
