# Literature and Nature-format audit

Checked 27 August 2026 against primary publications, DOI metadata and the
current Nature author guidance.

## Scientific framing

- **ARROW reference:** Pereyaslavets *et al.*, “Accurate determination of
  solvation free energies of neutral organic compounds from first principles,”
  *Nature Communications* (2022), DOI
  [10.1038/s41467-022-28041-0](https://doi.org/10.1038/s41467-022-28041-0).
  This is the source for the 85-solute comparison and the approximate 0.20
  kcal/mol PIMD8 result.
- **ARROW-NN:** Illarionov *et al.*, *JACS* (2023), DOI
  [10.1021/jacs.3c07628](https://doi.org/10.1021/jacs.3c07628). Its reported
  water-dimer interaction accuracy is not presented as a hydration-free-energy
  result.
- **NQE correction:** Kurnikov *et al.*, *JCTC* (2024), DOI
  [10.1021/acs.jctc.3c00921](https://doi.org/10.1021/acs.jctc.3c00921).
- **FreeSolv:** Mobley and Guthrie (2014), DOI
  [10.1007/s10822-014-9747-x](https://doi.org/10.1007/s10822-014-9747-x), and
  Matos *et al.* (2017), DOI
  [10.1021/acs.jced.7b00104](https://doi.org/10.1021/acs.jced.7b00104).
- **Continuum response:** SMD, DOI
  [10.1021/jp810292n](https://doi.org/10.1021/jp810292n); COSMO-RS, DOI
  [10.1021/j100007a062](https://doi.org/10.1021/j100007a062).
- **Solvation transfer:** CombiSolv, DOI
  [10.1016/j.cej.2021.129307](https://doi.org/10.1016/j.cej.2021.129307), and
  SoluteML, DOI
  [10.1021/acs.jcim.1c01103](https://doi.org/10.1021/acs.jcim.1c01103).

The manuscript makes no cross-protocol state-of-the-art claim. It distinguishes
the ARROW 85-solute reference set from a community benchmark and describes the
80/85 FreeSolv identity overlap only as an identity result, not label provenance.

## Closest prior art and novelty boundary

A second, focused audit on 28 August 2026 examined per-query solvent-conditioned
representations, surrogate-predicted physical descriptors and two-stage hydration
models. The closest architectural precedent is Jia *et al.*, *Chemistry--Methods*
(2026), DOI
[10.1002/cmtd.202500150](https://doi.org/10.1002/cmtd.202500150), which predicts 13
physically inspired solute descriptors from a molecular graph before hydration
regression. Stuyver and Coley's ml-QM-GNN previously established the same broad
two-stage pattern for reactivity, DOI
[10.1063/5.0079574](https://doi.org/10.1063/5.0079574).

The revised manuscript therefore makes no priority claim for surrogate-predicted
descriptors. It positions SolvAI's contribution as the controlled demonstration that
responses learned from calculated, empirical and corrected solvation sources provide
reusable molecular supervision. The complete map,
including ML-PCM, ImPerHam, 3D-RISM representations and systematic studies of when
predicted physical descriptors help, is recorded in
`reviews/NOVELTY_POSITIONING_AUDIT.md`.

## Nature Communications submission format

Current Nature Portfolio guidance was checked on 27 August 2026 and the manuscript
was subsequently rebuilt for *Nature Communications*. The release uses an explicit
abstract, numbered references, line numbering, figures integrated with the review
manuscript, a separate Methods section and a single Supplementary Information PDF
containing its supporting figures, tables and legends. A conservative standard LaTeX class is used because no
journal-specific class is required for initial submission.

The journal-specific declarations still requiring author confirmation are
isolated in `submission/AUTHOR_TODOS.md`.
