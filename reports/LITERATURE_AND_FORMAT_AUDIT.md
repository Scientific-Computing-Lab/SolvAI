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

## Nature submission format

The current [Nature formatting guide](https://www.nature.com/nature/for-authors/formatting-guide)
was checked on 27 August 2026. The release follows its practical submission
recommendations: a summary paragraph under approximately 200 words, numbered
references, line numbering, figures integrated with the text for review, a
separate Methods section and separate Supplementary Information. The source uses
a conservative standard LaTeX article class because no required public Nature
class is needed for initial submission.

The journal-specific declarations still requiring author confirmation are
isolated in `submission/AUTHOR_TODOS.md`.
