# SolvAI Figure 1 generation brief

## Scientific claim

Solvent-response information can be computed or measured once across external
datasets, distilled into structure-predictable molecular coordinates, and reused to
predict hydration free energy for a new molecule without running MD, PIMD, continuum,
alchemical, or probe calculations at inference time.

## Required full-figure narrative

Design a single integrated, full-width scientific Figure 1 with six clearly lettered
panels (A–F). It should read from upper left to lower right, but feel like one coherent
visual story rather than six unrelated boxes. The scientific hierarchy is:

1. Repeated physical simulation for every new solute is accurate but expensive.
2. External calculated, empirical, residual-corrected, and conformational sources
   expose complementary aspects of solvent response.
3. Six teacher families are distilled into 15 explicit response coordinates that can
   be predicted from molecular structure.
4. A separate endpoint model combines those 15 predicted priors with a 2,048-bit
   Morgan fingerprint and 217 RDKit descriptors, supervised by experimental hydration
   free energies.
5. Deployment is a simple learned path: SMILES to predicted response vector to
   hydration free energy. There is no MD, PIMD, or probe calculation at inference.
6. The matched experiment supports the claim: structure-only MAE 0.303; SolvAI MAE
   0.202 kcal mol^-1; reconstructed ARROW/PIMD8 reference MAE 0.205. PIMD8 is a
   separate accuracy reference, never a teacher or input to SolvAI. Shuffling the
   molecule-to-response alignment yields MAE 0.306 and abolishes the gain.

## Panel content

### A — The repeated-computation bottleneck

Show one chemically plausible small neutral organic solute moving from gas into a
water environment. Contrast a costly per-molecule trajectory/sampling route with the
goal of reusable supervision. Use a restrained visual cue for time/cost, not cartoon
money icons. Do not imply that SolvAI simulates water at inference.

### B — Complementary response teachers

Arrange six source families around a shared conceptual solvent-response space:

- CombiSolv-QM / COSMOtherm water response: 1 coordinate.
- Abraham empirical axes: 5 coordinates.
- OpenFF explicit-water response with residual correction: included within 2 corrected
  response coordinates together with GBn2.
- GBn2 implicit-water response with residual correction.
- MolSolv SMD(water): 1 coordinate.
- ConfSolv conformer reweighting and water-response distribution summaries: 6
  coordinates.

Do not depict these as equivalent-fidelity simulations and do not place their source
counts on a common quantitative axis. They deliberately span different approximations
and scales.

### C — Distillation into reusable coordinates

Show one neutral organic query molecule splitting into six independent frozen
structure-to-response surrogate families. The identical molecule must feed every
branch; do not use varied molecules, because molecule alignment is the scientific
point. Distinguish the actual surrogate types: CombiSolv-QM and SMD(water) use
CHEMELEON/D-MPNN; Abraham, OpenFF and GBn2 use ExtraTrees; ConfSolv uses LightGBM.
Each family produces only its own response coordinate group: CombiSolv-QM 1, Abraham
5, OpenFF 1, GBn2 1, SMD(water) 1 and ConfSolv 6. Concatenate these six outputs into
one compact, molecule-aligned 15-dimensional vector in the exact grouping
`1 | 5 | 1 | 1 | 1 | 6 = 15`. Never depict each surrogate as returning the entire
15-vector, and never route every family through one generic neural-network symbol.
State: "predicted per molecule by frozen structure→response surrogates · source
calculations are not run at inference". Emphasize that external benchmark-equivalent
molecules were excluded from supervised response teachers.

### D — Experimentally supervised endpoint

Show the predicted 15-response vector joining ordinary molecular structure features:
a 2,048-bit Morgan fingerprint plus 217 RDKit 2D descriptors. These enter an endpoint
ensemble trained against experimental hydration free energies. Keep the two training
stages visually distinct: response teachers are not trained on the ARROW/PIMD8 target,
and PIMD-derived features are absent from the final stack.

### E — Simulation-free deployment

Make this the visual focal point. A new SMILES string passes through the frozen
structure-to-response models and endpoint model to produce Delta G_hyd. Place a clear
small statement: "No MD · No PIMD · No probe calculation". Do not show water
trajectories inside the deployment path.

### F — Decisive matched evidence

Use a compact, honest quantitative inset rather than a decorative dashboard. Show:

- Matched structure-only: MAE 0.303 kcal mol^-1.
- SolvAI: MAE 0.202 kcal mol^-1.
- ARROW/PIMD8: MAE 0.205 kcal mol^-1, visibly separated and labeled "accuracy
  reference; not a teacher".
- Shuffled response priors: MAE 0.306 kcal mol^-1, showing that alignment is necessary.

Do not claim that SolvAI is statistically superior to PIMD8. The intended headline is
"PIMD8-level accuracy on this reference chemistry".

## Visual direction

Aim for the editorial clarity of a premium Nature or Science methods overview while
remaining original. Use a white or warm off-white background, abundant whitespace,
black typography, thin consistent connectors, and one accessible scientific palette.
Suggested semantic colors: deep navy for molecular structure, teal for learned
response coordinates, amber for physical source information, magenta only for the
separate PIMD8 reference, and neutral gray for excluded per-query simulation. Use
saturation only at the focal deployment/output path. Avoid gradients, glass effects,
glowing neural-network brains, generic AI circuitry, decorative laboratory icons,
3D text, drop shadows, and photorealistic fantasy chemistry. Molecular depictions must
remain modest and chemically plausible. Every arrow must have exactly one unambiguous
meaning. Keep all text horizontal and concise.

## Accuracy vetoes

- Never connect PIMD8 to training, response teachers, the 15 priors, or deployment.
- Never state or imply that the 15 responses are measured or simulated for each new
  query.
- Never depict the response priors as ground-truth decompositions of hydration free
  energy.
- Never imply that all six teacher families share a common fidelity or data scale.
- Never invent performance values, molecules, datasets, or response coordinates.
- Preserve the four reported MAEs exactly.
