# Michael--Boris collaborator follow-up audit

**Date:** 28 August 2026  
**Scope:** evidence and submission decisions only; the frozen models, folds and
manuscript were not changed.

## Executive assessment

| Item | Status | Answer in one line | Recommendation |
|---|:---:|---|---|
| 1. Identity, coverage and applicability | **A** | The final endpoint and every teacher are free of ARROW identities under exact and standardized-equivalence checks; teacher-dataset coverage is not an inference lookup requirement. | Supplement / Methods |
| 2. Which response sources matter? | **B** | The empirical/residual block is strongest alone; the narrow multi-source stack accounts for most of the gain, while SMD and ConfSolv add small, unresolved conditional increments. | Main + Supplement |
| 2a. What are the 15 response values? | **A** | They are 15 named scalar targets—empirical solute coordinates and calculated, corrected and conformational solvent responses—not opaque latent embeddings. | Collaborator clarification / Supplementary Table 1 |
| 3. Physically meaningful baselines | **B** | The causal baseline is 0.303 kcal mol\(^{-1}\); classical ARROW (0.785) supplies a legitimate physical scale, but the proposed historical 19-parameter model is not reproducible from the materials currently available. | Main; no approximate baseline |
| 4. Hard cases | **B** | The largest errors include methane, tert-butanol, two amides and several substituted/polycyclic aromatics; most large SolvAI misses are not shared by PIMD8. | Supplement + compact Discussion note |
| 5. \(\lambda\)/high-fidelity interpretation | **B** | Large structure-to-response error is real, but the learned \(dH/d\lambda\) and component curves are also protocol-conditioned rather than intrinsic molecular observables. | Discussion / Supplement |
| 6. Truly blind challenge | **C** | A valid challenge is possible only if an independent custodian freezes the molecules before predictions; it has not been run. Aspartic acid is unsuitable for the neutral-solute claim. | Optional prospective Supplement; do not delay submission |
| 7. Octanol transfer | **D** | Clean direct octanol data are too small for the proposed design, while the larger set is largely partition-derived; most SolvAI priors are water-specific. | Future work / separate study |
| 8. Closest literature | **A** | The current Introduction already credits the closest two-stage hydration, per-query physical representation, computed-transfer and surrogate-descriptor precedents. | No action |

Here, **A** means already answered convincingly, **B** means a new analysis of
frozen outputs, **C** means a clean prospective run would be required, and **D**
means outside the present manuscript. The only new executable analysis was
descriptive:

```bash
uv run python scripts/analyze_collaborator_followup.py
```

It fits no model and writes machine-readable summaries to
`results/followup_audit/`. All endpoint values below use the conservative
standardized-exclusion teachers and frozen confirmatory predictions.

## Quick collaborator Q&A

### Does a new molecule need to exist in all six teacher datasets?

**No.** Those datasets were used to train frozen structure-to-response surrogates;
at inference, each surrogate predicts its response coordinate directly from a new
SMILES rather than looking the molecule up in its source table. This is computational
coverage, not proof of unrestricted validity: the tested domain is neutral,
predominantly small organic hydration chemistry, not ions, salts, metals,
biomacromolecules or chemically remote structures.

### Does SolvAI currently claim applicability to charged or ionic solutes?

**No.** The paper validates hydration of **neutral organic solutes**: all 85 aqueous
ARROW benchmark entries have formal charge zero, and the endpoint is trained with
1,280 benchmark-disjoint neutral experimental hydration labels plus neutral ARROW
outer-training rows. The manuscript and model card explicitly exclude validation for
ions, salts, metals, proteins and broad chemical extrapolation. A parseable charged
SMILES may still produce a numerical output, but that is software behaviour—not an
evidence-backed prediction claim.

Accordingly, Li\(^+\), Cl\(^-\), glutamate\(^{2-}\), charged or zwitterionic
aspartate, guanidinium and the usual aqueous zwitterionic form of GABA are
**out-of-scope prospective stress tests**. Dipeptides are likewise outside the
validated domain because of their peptide-like size, flexibility and commonly
zwitterionic microstates; even a neutral or capped representation would test transfer
limits rather than established applicability. Neutral, single-component substituted
aromatics are closer to the current domain and can be useful challenging
generalization tests, but an unseen example remains prospective and the frozen audit
already shows that some substituted/polycyclic aromatics are difficult. None of
these tests should be folded into the present neutral-solute accuracy claim or used
post hoc as a favourable demonstration. The controlling evidence is the
[manuscript Methods](https://github.com/Scientific-Computing-Lab/SolvAI/blob/main/paper/main.tex),
[Supplementary Methods](https://github.com/Scientific-Computing-Lab/SolvAI/blob/main/paper/supplementary/supplementary.tex),
and [model card](https://github.com/Scientific-Computing-Lab/SolvAI/blob/main/models/final/MODEL_CARD.md).

### Does the advantage survive genuinely different chemistry?

**Yes as a relative advantage, but not at the same absolute accuracy.** When chemical
separation was applied to *all* endpoint-labelled training molecules, family holdout
gave **0.468 versus 1.262**, scaffold holdout **0.376 versus 0.728**, and a Butina
cluster split **0.216 versus 0.349** kcal mol\(^{-1}\) for SolvAI versus the matched
structure-only model. After excluding every training molecule above Morgan
similarities 0.50, 0.60, 0.70 and 0.80 to each test molecule, the corresponding pairs
were **0.219/0.359, 0.211/0.344, 0.211/0.343 and 0.210/0.339**. The priors therefore
retain value under genuine separation, while the family/scaffold results also show
that broad chemical extrapolation remains harder than interpolation.

### Which source contributes most—is one source responsible?

The strongest block *in isolation* is the empirical/residual-corrected set (Abraham
axes plus corrected OpenFF/GBn2): **0.243 versus 0.303** kcal mol\(^{-1}\) for
structure only. Computation-derived core and SMD-only blocks each give about **0.253**;
ConfSolv alone is nearly neutral at **0.299**. In combination, however, the narrow
multi-source stack reaches **0.212**, adding SMD gives **0.207**, and the full stack
gives **0.202**. The small conditional SMD and ConfSolv increments are not separately
resolved by their paired intervals, so the evidence supports complementarity across
an aligned response layer—not a claim that one teacher causes the full result.

### Why did the λ/alchemical-response experiment fail?

Two limitations coincide. First, the structure-to-response heads were not accurate
enough: their MAEs were **1.275–5.199 kcal mol\(^{-1}\)**, and adding the predicted
PIMD2 responses worsened the matched endpoint from **0.19592 to 0.20136**. Second,
intermediate \(dH/d\lambda\) values are conditioned on the chosen coupling schedule,
soft-core form and parameterization; unlike the integrated free-energy endpoint,
their detailed curve is not a unique intrinsic molecular target. Sparse supervision
and protocol-dependent targets therefore both limit transfer.

### Where does SolvAI fail?

The largest frozen errors are methane (**1.102**), tert-butanol (**1.052**),
N-methylacetamide (**1.040**), m-cresol (**0.810**), tetracene (**0.754**),
dimethoxyethane (**0.723**), pyrene (**0.717**) and dimethylacetamide (**0.700**)
kcal mol\(^{-1}\). Amides and larger/conjugated aromatics are represented among the
hard cases, but they are not the whole pattern: methane and tert-butanol are the two
largest misses, and amides/aromatics improve on average relative to structure only.
Of nine SolvAI errors at least 0.5 kcal mol\(^{-1}\), only m-cresol is also a
PIMD8 error above that threshold, so most are SolvAI-specific rather than universally
hard reference molecules.

### Are there already methods doing something similar?

**Yes, for important parts of the idea.** Prior work includes computed-to-experimental
transfer, two-stage hydration models, per-query physical representations such as
ML-PCM, ImPerHam and 3D-RISM, and surrogate-predicted QM descriptors. SolvAI is not
claimed as the first two-stage or predicted-descriptor model. Its supported
distinction is the controlled use of a benchmark-disjoint, multi-source
solvent-response layer whose value survives matched removal, molecule-wise shuffling,
repeated partitions, global chemical separation and zero-ARROW-label transfer, while
the originating calculations are absent at inference. PIMD8 is the accuracy
comparator, not a retained teacher.

## 1. Identity, coverage and applicability — A

### Collaborator question

Could common ARROW molecules have appeared under another identity in an external
endpoint or teacher dataset, and does incomplete source coverage prevent SolvAI from
predicting some new molecules?

### What the question tests

This tests both label leakage and a separate deployment concern. Leakage asks
whether a held-out ARROW molecule supplied supervised information to either learning
stage. Coverage asks whether inference requires finding the query molecule in every
teacher's original source table.

### Existing evidence

The final audit applies canonical isomeric SMILES, full InChIKey and the first
InChIKey block, followed by fragment-parent, uncharged-parent and canonical-tautomer
standardization. The 1,280-label endpoint pool contains **zero** ARROW matches under
all of those definitions. This is important because **80 of the 85** ARROW
connectivities occur in FreeSolv; those identities were excluded rather than treated
as external endpoint labels.

The original connectivity filter removed the following source records before
teacher fitting:

| Teacher source | After exact filter, before standardized refit | Exact ARROW overlaps removed |
|---|---:|---:|
| CombiSolv-QM water | 3,963 | 25 |
| SoluteML Abraham | 8,098 | 84 |
| OpenFF explicit-water ASFE | 520 | 83 |
| GBn2 implicit solvent | 550 | 0 |
| MolSolv SMD(water) | 350,391 | 82 |
| ConfSolv H2O | 39,878 | 13 |

The expanded equivalence audit then found records missed by exact connectivity:
**2 CombiSolv-QM**, **32 MolSolv** and **22 ConfSolv** rows. Each was removed, the
three teachers were refitted while preserving all remaining source split membership,
and the complete endpoint, repeat, shuffled-prior and chemical-separation analyses
were rerun. Verification passed. No standardized-equivalent endpoint row remained.
The final fixed teacher splits contain 3,959 CombiSolv-QM rows and 350,359 MolSolv
rows; the selected complete-target ConfSolv fitting table contains 17,829 rows.
Morgan-similarity-1 collisions between different connectivities were retained as
fingerprint collisions and explicitly listed; they were not relabelled as identities
or silently removed.

Primary records are
`audits/leakage_audit.{json,csv,md}`,
`audits/confirmatory/standardized_exclusion_records.csv`,
`audits/confirmatory/standardized_exclusion_refit_verification.json`, and
`audits/confirmatory/chemical_identity_matches.csv`.

### Coverage at inference

Source coverage is used to **fit** a response surrogate; it is not a database lookup
at deployment. For a new parseable SMILES, SolvAI computes Morgan/RDKit structure
features, runs each frozen structure-to-response model and obtains all 15 response
coordinates. The query therefore need not occur in MolSolv, ConfSolv, CombiSolv or
another teacher table. The final artifact contains no source-coverage gate.

This computational ability is not evidence of unrestricted chemical validity. The
validated domain is neutral, predominantly small organic hydration chemistry similar
to the source and reference domains. Ions, salts, metals, biomacromolecules, ambiguous
aqueous microstates and substantially larger or chemically remote molecules are not
validated merely because RDKit can parse their SMILES.

### Exact result and manuscript implication

The common-molecule concern is answered: supervised ARROW identity overlap is zero
under the final conservative policy. The coverage concern should be explained with
one sentence distinguishing teacher training coverage from inference availability,
followed by the domain limitation.

**Recommendation:** Supplement / Methods. No new run and no main-text expansion are
needed.

## 2. Which response sources matter? — B

### Collaborator question

Which physical-response sources most influence the endpoint, and is one teacher
responsible for the result?

### What the question tests

This asks four different questions that should not be conflated: how accurately a
surrogate predicts its own source target; what a block contributes by itself; what it
adds conditional on other priors; and whether apparently useful sources are
redundant.

### New analysis of frozen outputs

The frozen matched block results were consolidated and paired molecule bootstrap
intervals were recomputed without fitting or selecting a model. The blocks are
predeclared scientific ablations, not an additive variance decomposition.

| Endpoint input beyond structure | MAE | Change from 0.303 structure baseline | 95% paired interval | Interpretation |
|---|---:|---:|---:|---|
| Empirical/residual-corrected: Abraham + corrected OpenFF/GBn2 | 0.243 | -0.061 | [-0.108, -0.022] | Strongest isolated block |
| Computation-derived core: CombiSolv-QM + raw OpenFF/GBn2 | 0.253 | -0.051 | [-0.144, 0.029] | Lower point estimate; interval crosses zero |
| SMD(water) alone | 0.253 | -0.050 | [-0.128, 0.006] | Lower point estimate; interval crosses zero |
| ConfSolv six-response block alone | 0.299 | -0.004 | [-0.014, 0.004] | Neutral alone |
| Narrow cumulative stack | 0.212 | -0.092 | [-0.183, -0.019] | Most of the full gain is already present |
| Narrow stack + SMD | 0.207 | -0.097 | [-0.209, -0.016] | Best point improves by 0.005 over narrow |
| Full 15-prior SolvAI | 0.202 | -0.101 | [-0.215, -0.020] | Best predeclared point estimate |

The incremental SMD change relative to the narrow stack is -0.005 kcal mol\(^{-1}\)
with a paired interval of [-0.032, 0.015]. Adding ConfSolv after that changes MAE by
-0.004, interval [-0.010, 0.002]. Thus the endpoint data do not resolve either small
conditional increment separately. They support a large effect for the response layer
as a whole and for the narrow block, not a causal ranking of every teacher.

The destructive control is decisive: aligned full priors give **0.202** MAE, whereas
five molecule-wise shuffled versions average **0.307**, close to structure only.
Across five frozen partitions, structure only is **0.313 ± 0.004** and full SolvAI is
**0.207 ± 0.004** kcal mol\(^{-1}\), with full SolvAI better on every partition.

### Teacher fidelity is a different quantity

| Source | Own-target validation |
|---|---|
| CombiSolv-QM | 0.739 kcal mol\(^{-1}\) MAE, fixed source test (n=395) |
| MolSolv SMD(water) | 0.712 kcal mol\(^{-1}\), fixed source test (n=17,520) |
| Abraham E/S/A/B/L | 0.035–0.240 on their distinct Abraham scales, source OOF |
| OpenFF | 0.617 for raw \(\Delta G\), 0.531 for experimental residual, source OOF |
| GBn2 | 0.669 for raw \(\Delta G\), 0.626 for experimental residual, source OOF |
| ConfSolv selected targets | 0.248–0.688 kcal mol\(^{-1}\), fixed source validation (n=1,785) |

These numbers differ in units, scale, target variance and validation protocol. They
cannot be ranked as if lower teacher MAE implied greater endpoint value. ConfSolv is
the clearest counterexample: several targets are predicted accurately on their own
scale, yet its isolated endpoint block is neutral.

### Scientific interpretation

The most transferable information is not simply the highest-fidelity calculation.
The strongest isolated signal combines empirical solute-response axes with
experimental residual corrections to approximate water calculations. Most of the
endpoint gain emerges from a compact, mixed response vocabulary; SMD and ConfSolv
then provide small, statistically unresolved conditional refinements. This supports
**relevance + structure-learnability + complementarity**, not a claim that one source
dominates.

Machine-readable results are
`results/followup_audit/source_block_summary.csv` and
`results/followup_audit/teacher_fidelity_summary.csv`.

**Recommendation:** show the matched block result and shuffled control in the main
paper; place full teacher-fidelity detail in Supplementary Information. Do not show a
raw tree-importance ranking as causal evidence.

## 2a. What exactly are the 15 learned response values? — A

### Collaborator question

What physical or physicochemical quantity does each response teacher predict?

### Answer

The 15 values are **not opaque latent embeddings**. Each is a predefined scalar target
with a named interpretation, learned from molecular structure and predicted for a new
SMILES by a frozen surrogate:

| Intuitive group | Values | Meaning |
|---|---|---|
| Polarity and hydrogen bonding | **5 Abraham axes:** E, S, A, B and L | Excess molar refraction; dipolarity/polarizability; hydrogen-bond acidity; hydrogen-bond basicity; and hexadecane--air partition response. These are empirical physicochemical coordinates, not pure physical observables. |
| Calculated bulk-water solvation | **2:** CombiSolv-QM/COSMOtherm water and MolSolv SMD(water) | Predicted water solvation free energies from two distinct calculated solvation formalisms (kcal mol\(^{-1}\)). |
| Explicit- and implicit-water response | **2:** corrected OpenFF and corrected GBn2 | Respectively, an explicit-water alchemical and an implicit-solvent hydration prediction, each **plus a separately learned experimental residual** (kcal mol\(^{-1}\)); these corrected coordinates are therefore not pure simulation observables. |
| Conformational water response | **6 ConfSolv summaries** | Gas conformer correction; solution conformer correction; hydration conformer correction; conformer solvation-energy spread; and the mean and spread of the solvent-induced change in conformer relative energies (all kcal mol\(^{-1}\)). |

Thus the layer is best described as an interpretable mixture of **empirical,
calculated, residual-corrected and conformational response coordinates**. The exact
definitions are frozen in [Supplementary Table 1](https://github.com/Scientific-Computing-Lab/SolvAI/blob/main/paper/supplementary/tables/response_priors.tex),
and their assembly into the 15-value inference vector is explicit in
[`solv_ai/teachers.py`](https://github.com/Scientific-Computing-Lab/SolvAI/blob/main/solv_ai/teachers.py).
No PIMD, NQE or \(\mathrm{d}H/\mathrm{d}\lambda\) value is among the retained 15.

**Recommendation:** collaborator clarification; Supplementary Table 1 already carries
the definitive per-coordinate specification. No new analysis or manuscript change is
required for this audit.

## 3. Accuracy on a physical baseline scale — B

### Collaborator question

Can the reader see the physical error scale before SolvAI, including the older
approximately 19-parameter solvent-accessible/generalized-Born-style model discussed
in the meeting?

### Existing evidence

Two comparisons answer distinct questions:

- **Causal model comparison:** with identical 1,280 external experimental labels,
  fold-local ARROW labels, weights, Morgan/RDKit representation, ExtraTrees settings,
  folds and seeds, removing only the priors changes MAE from **0.202 to 0.303**.
- **Physical context on the same 85 solutes:** reconstructed classical ARROW is
  **0.78465** MAE and ARROW/PIMD8 is **0.20484** kcal mol\(^{-1}\).

The frozen source for these values is `results/paper_metrics.json`.

The 0.785 value is useful context for the physical accuracy scale, but it is not a
matched ML ablation. Conversely, 0.303 is the correct causal baseline and should not
be replaced by a deliberately weak historical model.

### Audit of the proposed older model

The complete release, exploratory workspace and recovered Freecurve repositories
were searched for an identifiable 19-parameter implementation, parameter table,
atom-typing rules or frozen predictions. None was found. Available materials include
generic GBSA implementations, a 562-molecule published GBSA-AMBER table and an
exploratory continuum-polarization/SASA teacher. They are not the proposed historical
model. The latter was tested only as **surrogate-predicted features** in an earlier,
different endpoint pipeline (0.2523 versus a 0.2504 matched direct model); it is not a
standalone physical prediction on the current 85-solute evaluation and cannot supply
the requested comparator.

An apples-to-apples reproduction presently lacks the exact citation/model version,
19 fitted parameters, atom typing, charge and conformer protocol, standard-state
convention and output sign convention. Choosing substitutes after seeing the ARROW
targets would manufacture a baseline.

### Exact result and manuscript implication

No new historical MAE is reportable. If Michael supplies the exact frozen
implementation or complete specification, it can be run once prospectively on all
compatible neutral ARROW molecules, with coverage and exclusions reported. Until
then, the manuscript should use 0.303 as its primary matched baseline and 0.785
classical ARROW as physical context.

**Recommendation:** Main for the 0.303 and 0.785 scales; no approximate 19-parameter
bar. Adding a future result is non-blocking and would belong in Supplementary
Information unless exact comparability is unusually strong.

## 4. Hard cases and failure modes — B

### Collaborator question

Are failures concentrated in substituted aromatics, larger/conjugated systems,
amides or another chemically interpretable class, and are the same molecules hard for
PIMD8?

### New analysis of frozen predictions

The table uses only the frozen standardized-exclusion OOF predictions. It is
post-hoc diagnosis, not confirmatory hypothesis testing.

| Rank | Molecule | SMILES | Family | SolvAI error | Structure-only error | PIMD8 error |
|---:|---|---|---|---:|---:|---:|
| 1 | Methane | `C` | Alkane | 1.102 | 0.691 | 0.091 |
| 2 | tert-Butanol | `CC(C)(C)O` | Alcohol | 1.052 | 0.058 | 0.350 |
| 3 | N-Methylacetamide | `CNC(C)=O` | Amide | 1.040 | 1.292 | 0.494 |
| 4 | m-Cresol | `Cc1cccc(O)c1` | Aromatic | 0.810 | 0.879 | 1.029 |
| 5 | Tetracene | `c1ccc2cc3cc4ccccc4cc3cc2c1` | Aromatic | 0.754 | 0.771 | 0.440 |
| 6 | Dimethoxyethane | `COCCOC` | Ether | 0.723 | 1.148 | 0.197 |
| 7 | Pyrene | `c1cc2ccc3cccc4ccc(c1)c2c34` | Aromatic | 0.717 | 0.817 | 0.011 |
| 8 | Dimethylacetamide | `CC(=O)N(C)C` | Amide | 0.700 | 1.447 | 0.009 |
| 9 | Acetic acid | `CC(=O)O` | Acid | 0.551 | 0.206 | 0.342 |
| 10 | Phenol | `Oc1ccccc1` | Aromatic | 0.488 | 0.854 | 0.181 |

Nine molecules have SolvAI absolute error at least 0.5 kcal mol\(^{-1}\). Only
m-cresol also has PIMD8 error at least 0.5; the other eight are predominantly
SolvAI-specific misses under this descriptive threshold.

The meeting hypotheses receive qualified support:

- **Amides:** n=4, SolvAI MAE 0.483 versus 0.708 structure only and 0.184 PIMD8.
  The priors help the class overall, but two amides remain among the eight worst
  individual predictions.
- **Aromatic structures:** n=17, 0.265 versus 0.426 structure only and 0.317 PIMD8.
  For the post-hoc subset with an aromatic ring and at least eight heavy atoms
  (n=10), the values are 0.334, 0.430 and 0.388, respectively. Large/conjugated
  aromatics are enriched among hard cases, but the sample is too small for a ranked
  population claim.
- **Other failures matter:** methane and tert-butanol are the two largest errors,
  and the priors worsen the average for alkanes (0.234 versus 0.201) and alcohols
  (0.196 versus 0.109). The failure mode is not simply molecular size.

Across all 85 molecules, SolvAI improves absolute error over the matched structure
model for 53. The broad gain is therefore not produced by every molecule, and the
largest misses are not generally inherited from the PIMD8 comparator.

Full records, including scaffolds, descriptors, signed residuals and five-repeat
stability, are in `results/followup_audit/hard_case_audit.csv`; grouped summaries are
in `hard_case_group_summary.csv`.

**Recommendation:** Supplement for the molecule-level table/plot; one compact
Discussion sentence identifying amides and aromatic systems while noting methane and
tert-butanol. Do not promote the post-hoc group ordering to a formal generalization
claim.

## 5. Failed \(\lambda\) and high-fidelity response — B

### Collaborator question

Did direct high-fidelity response distillation fail only because the surrogate was
inaccurate, or because the intermediate target is itself tied to a chosen alchemical
protocol?

### What was actually learned

The experiment used 5-ps, two-bead PIMD observations at \(\lambda=0.1,0.5,0.9\):
total \(dH/d\lambda\), Coulomb and van der Waals contributions. There were 72, 74 and
73 successful molecules at the three states, with 72 complete curves. The exported
polarization component was identically zero. The ARBALEST configuration used a
specified 15-state coupling schedule, scale-factor power 2.0, soft-core power 1.0,
radius 1.5 and cutoff 6.0. These settings are visible in the archived `conf.xml`
files; the summary is in `reports/LAMBDA_RESPONSE_EXPERIMENT.md`.

### Exact result

Response-head MAEs span **1.275–5.199 kcal mol\(^{-1}\)**. The matched endpoint changes
from **0.19592** to **0.20136** with predicted PIMD2 response, to **0.21474** with the
classical/NQE/PIMD hierarchy, and to **0.21901** with both. Direct integration of the
predicted three-point curve followed by fold-local affine calibration gives **1.51356**
kcal mol\(^{-1}\).

### Technical interpretation

Path dependence applies. \(dH/d\lambda\) is the derivative of a selected interpolating
Hamiltonian. Its shape—and especially its electrostatic/van der Waals decomposition—
depends on the coupling map, soft-core treatment, parameterization and component
convention. A converged free-energy difference between fixed physical endpoints is a
state function; a sparse short-run three-point curve and its components are not unique
intrinsic molecular observables.

The existing explanation is therefore correct but incomplete. The measured
structure-to-response errors are too large for the 0.2-kcal mol\(^{-1}\) endpoint, but
the difficulty reflects both sparse/noisy supervision and a protocol-conditioned
target. This experiment does **not** show that all alchemical-response distillation is
unhelpful.

Proposed wording for later manuscript consideration:

> The protocol-specific short-run responses were learned with 1.27–5.20 kcal
> mol\(^{-1}\) errors from only 72–74 labels and did not improve the endpoint. Their
> detailed shapes depend on the alchemical path, force field and soft-core convention;
> the result therefore bounds this sparse response representation rather than
> alchemical-response distillation in general.

**Recommendation:** Discussion / Supplement. This is a scientifically material
interpretive correction, but no new experiment is required.

## 6. Truly blind molecule challenge — C

### Collaborator question

Can an independently selected difficult molecule provide a genuinely blind test?

### Feasibility and scientific value

Yes, but only prospectively. Choosing a molecule now from known errors would be
anecdotal selection. A single molecule would be a demonstration, not a statistically
meaningful validation; a small panel is much stronger.

Aspartic acid should not be used automatically. In water near neutral pH it is
predominantly zwitterionic/anionic, with multiple protonation microstates. A formally
neutral covalent SMILES does not make its experimental hydration free energy a clean
neutral-nonelectrolyte endpoint. It is outside the present validated scope unless a
specific microstate, pH-independent thermodynamic convention and matching experimental
quantity are established in advance.

### Exact prospective protocol

1. Boris or another independent custodian freezes **20–30** candidate SMILES and
   experimental references before receiving any SolvAI prediction. A single showcase
   molecule may accompany this panel but cannot replace it.
2. Eligibility is applied outcome-blind: one neutral fragment; no metal; no ambiguous
   dominant aqueous protomer; neutral small-organic elements and size within a
   prespecified reference/training envelope; a compatible 298-K gas-to-water standard
   state; and a stated experimental uncertainty.
3. Before inference, audit each candidate against the 1,280 endpoint rows and every
   teacher source using full identity, connectivity, fragment parent, uncharged parent
   and canonical tautomer. Predeclare how ineligible candidates are replaced.
4. Freeze candidate identities, source citations, endpoint values under custodian
   holdback, the SolvAI artifact hash and the analysis script. Record scaffold and
   nearest-neighbour similarity without using them to select on performance.
5. Run one batch prediction, unblind all ground truth and report every molecule,
   including exclusions and failures.

No challenge was run because the necessary independent preselection has not occurred.

**Recommendation:** optional prospective Supplement if it can be conducted cleanly
and quickly; otherwise future work. Submission should not wait for a single-molecule
demonstration.

## 7. Different solvent, especially octanol — D

### Collaborator question

Would transfer to octanol demonstrate that the response-learning framework extends
beyond hydration?

### Existing endpoint data

The most relevant local source is
[SoluteML](https://doi.org/10.1021/acs.jcim.1c01103) (Zenodo 5792296). For
1-octanol:

| SoluteML table | All-data unique solutes | Selected neutral unique solutes | Underlying entries |
|---|---:|---:|---:|
| dGsolvDB1 / dGsolvDB2 direct compiled values | 260 | 252 | 397 all / 389 selected |
| dGsolvDB3 expanded values | 4,304 | 4,261 | 8,980 all / 8,911 selected |
| logPow table | 3,061 | 3,019 | partition coefficients, not direct octanol \(\Delta G\) |

The apparently large dGsolvDB3 subset does not provide 4,304 independent direct
octanol solvation measurements. The source paper states that dGsolvDB3 augments
dGsolvDB2 by converting water–octanol and other solvent–water partition coefficients
to solvation free energies under a dry-solvent approximation; the octanol rows are
dominated by in-house, PHYSPROP and OCHEM logP/logPow sources. This derived provenance
creates correlated labels and a water-reference dependency that must be treated
explicitly. The clean direct subset is far below the proposed roughly 1,000 training
plus 100 held-out cases.

### Compatibility with SolvAI priors

Ten of the 15 deployed priors are explicitly water-specific: CombiSolv-QM water,
corrected OpenFF explicit water, corrected GBn2 water, SMD(water), and six ConfSolv
H2O summaries. The five Abraham solute axes are broadly solvent-transferable, but an
octanol endpoint requires solvent-specific coefficients or learning. Reusing the ten
water priors may be a valid cross-solvent representation experiment, but it would not
mean that SolvAI had distilled octanol response. A persuasive extension would require
curating an octanol endpoint with dependency-aware splitting and training
octanol-specific response teachers—or explicitly narrowing the claim to transfer from
water response.

### Go/no-go for this submission

**No-go.** Endpoint-only refitting would answer a different and weaker question;
building solvent-specific teachers and a clean independent test is a substantial new
study. It risks obscuring the controlled hydration result and should not delay the
current submission.

**Recommendation:** Future work / separate paper on solvent-conditioned response
distillation.

## 8. Literature and “does anyone else do this?” — A

### Collaborator question

Does the manuscript acknowledge methods that already combine predicted or computed
physical descriptors with molecular learning, and is the remaining contribution
clear?

### Existing evidence

The current Introduction now gives full credit to the closest lines:

- direct hydration learning and database-bias studies;
- computed-to-experimental transfer, including CombiSolv;
- single-calculator solvation surrogates;
- physically inspired hydration descriptors;
- per-query ML-PCM, 3D-RISM, ImPerHam and hybrid alchemical/ML representations;
- surrogate-predicted QM descriptors (ml-QM-GNN);
- controlled studies of when predicted QM descriptors help;
- surrogate-descriptor versus latent-representation studies; and
- the 2026 two-stage physically inspired hydration framework, the closest direct
  architectural precedent.

The text explicitly says that these works establish the architectural ingredients.
It then states the narrower unresolved scientific question that the confirmatory
controls actually answer: whether a molecule-aligned layer assembled from calculated,
empirical and corrected solvation-response sources contributes beyond an otherwise
identical experimentally supervised endpoint after source calculations are removed
from inference.

### Exact conclusion

The answer to “are there already other methods doing this?” is **yes for individual
ingredients and for two-stage predicted descriptors**, and the manuscript now says so.
SolvAI's supported contribution is the controlled multi-source solvent-response test,
not priority for two-stage learning. Its matched removal, molecule-response shuffling,
five partitions, global chemical separation and zero-ARROW transfer distinguish the
evidence. PIMD8 remains an accuracy comparator, not a retained teacher.

No missing comparator was found that materially overturns this positioning. The
closest-paper audit is preserved in `reviews/NOVELTY_POSITIONING_AUDIT.md` and
`reviews/NOVELTY_EDITORIAL_CHECK.md`.

**Recommendation:** No action. Reopening the literature section now would add volume,
not change the claim.

## Collaborator-facing response

We checked all eight questions against the frozen models and outputs, and ran only
descriptive analyses where the existing tables did not already give a direct answer.
The identity firewall is stronger than exact-SMILES exclusion: all ARROW-equivalent
endpoint records are absent, standardized tautomer/parent matches were removed from
three teacher sources, and those teachers and all confirmatory endpoints were refit.
Teacher-dataset coverage does not impose an inference lookup requirement—the frozen
surrogates produce every response prior from a new SMILES—but validation remains
limited to neutral small-organic chemistry.

The response gain is real but not attributable to one magic teacher. The strongest
isolated block is the empirical/residual-corrected response set (0.243 versus 0.303),
the compact narrow stack supplies most of the improvement (0.212), and the complete
15-prior model reaches 0.202; shuffling the priors returns performance to 0.307. SMD
and ConfSolv provide small conditional point improvements that are not separately
resolved. The worst SolvAI cases include methane, tert-butanol, amides and several
aromatics; only one of the nine errors above 0.5 kcal mol\(^{-1}\) is also above 0.5
for PIMD8.

The proposed old 19-parameter baseline cannot be reproduced without its exact model
specification, so we recommend showing the matched 0.303 control and classical ARROW
at 0.785 rather than inventing a substitute. The failed \(\lambda\) experiment should
be described more precisely: its errors are large, and its sparse component curves
are protocol-dependent rather than unique molecular ground truth. A blind challenge
is worthwhile only if Boris freezes a panel before prediction; aspartic acid is not a
clean neutral-solute test. Octanol is a separate project: the clean direct dataset is
only about 252 neutral solutes, the larger set is largely partition-derived, and ten
of the 15 current priors are water-specific. We therefore recommend no octanol
extension before submission and no further literature rewrite.
