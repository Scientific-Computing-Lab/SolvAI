# Cover letter

Dear Editors,

We submit “Structure-predicted solvent responses enable simulation-free hydration
free-energy prediction” for consideration as an Article in *Nature Communications*.

High-accuracy solvation calculations repeatedly pay the cost of resolving solvent
response for each molecule. Prior work has shown that physical descriptors can enrich
molecular learning, and that descriptors normally obtained from quantum calculations
can themselves be predicted from structure. SolvAI advances this line of research by
asking whether complementary solvent responses learned from calculated, empirical and
corrected sources can become reusable molecular supervision. It predicts those
responses from structure and exposes them to a separately supervised experimental
endpoint. A deployed prediction accepts a SMILES string and runs no molecular
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

We also prospectively qualified an external cohort before evaluating either endpoint.
All 220 retained molecules are disjoint from every experimental endpoint label;
97 are additionally absent from all six response-teacher source tables. On these two
cohorts the same matched contrast lowers MAE from 1.532 to 1.153 and from 2.138 to
1.536 kcal mol⁻¹, respectively, with both paired intervals excluding zero. The higher
absolute errors expose the current domain limit while showing that the response-layer
advantage is not confined to the ARROW molecules or direct teacher-source exposure.

Building on established two-stage architectures, the work provides a controlled
demonstration that aligned solvent responses carry complementary endpoint information,
together with an experimentally defined boundary: compact, learnable responses
transfer, whereas sparse high-fidelity alchemical responses currently do not because
their structure-to-response errors remain too large. This distinction makes the study
relevant beyond hydration prediction. It shows how costly physical-response data can
be converted into reusable molecular supervision, while establishing when that
strategy is likely to succeed.

The manuscript is accompanied by molecule-level predictions, the preregistration,
identity and similarity audits, all negative controls, frozen artifacts, source code
and a complete reproducibility package.

Sincerely,

Gal Oren and Michael Levitt
