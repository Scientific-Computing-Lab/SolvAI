# SolvAI Figure 1: scientific asset and design record

## Selected composition

Variant C, **balanced**, is the publication figure. Its editable source of truth is
`SolvAI_Figure1.penpot` (SHA-256
`92a6933f93eac6fd7dbad39bec2a23026f19363398061a7872d94c6fff7ee9c3`).
The file contains the design system and all three alternatives on separate pages.
All layout, type, rules, arrows, response-coordinate outlines and model glyphs are
native Penpot objects. Imported scientific objects remain editable SVG groups.

## Science-first asset decisions

- **Molecular structure:** N-methylacetamide (`CNC(C)=O`) rendered by RDKit.
  Indigo was tested independently; RDKit was retained because its depiction was
  clearer at journal size.
- **Calculated solvation:** a molecular cavity built deterministically from the
  RDKit coordinates and atomic van der Waals radii. It is an explanatory cavity
  schematic, not an electrostatic-potential calculation.
- **Water response:** a fixed Packmol configuration containing N-methylacetamide
  and 72 waters; the selected view is a vector projection of the nearest 18 waters.
  It is illustrative geometry, not an equilibrated trajectory or model input.
- **Conformational response:** three actual RDKit ETKDGv3/MMFF conformers of
  1,2-dimethoxyethane, generated with seed `0x5A17` and retained with coordinates.
- **Quantitative diagnostic:** variant B uses real N-methylacetamide PIMD2
  observations at lambda 0.1, 0.5 and 0.9. It is explicitly marked explored and
  not retained; it is absent from the selected figure.
- **External libraries:** SciDraw and Bioicons were searched before composition.
  Two CC0 Bioicons interaction glyphs are retained only in variant B. SciDraw's
  water-drop artwork was rejected because the Packmol-derived solvent geometry was
  scientifically more faithful. URLs, licenses and modifications are recorded in
  `fig1_assets/EXTERNAL_ASSET_MANIFEST.json`.

## Response-coordinate audit

The selected panel contains exactly the 15 scalar outputs defined in
`paper/supplementary/tables/response_priors.tex`:

1. COSMOtherm water response
2. Abraham E
3. Abraham S
4. Abraham A
5. Abraham B
6. Abraham L
7. corrected OpenFF explicit-water response
8. corrected GBn2 implicit-solvent response
9. SMD(water) response
10. gas conformer correction
11. solution conformer correction
12. hydration conformer correction
13. conformer solvation-energy spread
14. mean conformer solvent response
15. conformer response spread

These are named, interpretable scalar targets rather than a latent embedding. The
Abraham coordinates are empirical physicochemical axes; the corrected OpenFF and
GBn2 coordinates include learned experimental residuals.

## Architecture audit

The figure shows the released stack accurately:

- Stage 1: six frozen structure-to-response surrogate artifacts produce 15 values.
- Stage 2: 2,265 deterministic structure features (2,048 Morgan bits plus 217 RDKit
  descriptors) and the 15 predicted responses feed a three-member ExtraTrees
  endpoint ensemble; experimental hydration free energies supervise this endpoint.
- Deployment: a SMILES string drives both the deterministic feature path and the
  frozen response surrogates. The original source calculations do not run.
- PIMD2, PIMD4, PIMD8 and lambda response are not among the retained 15 values.
  ARROW/PIMD8 is shown only as the separate accuracy reference.

## Variants and selection

- **A — minimal:** fastest conceptual read, but deliberately suppresses much of the
  source chemistry.
- **B — molecular-science-forward:** strongest molecular/interaction detail, but
  gives the explored, non-retained lambda response more visual weight.
- **C — balanced:** selected because it preserves real molecular evidence, makes
  the four response families and 15 named coordinates explicit, and gives the
  simplest uninterrupted path to deployment.

All variants were rendered and inspected at 180-mm width. The contact sheet is
`previews/fig1_penpot_variants_contact.png`.

## Reproduction and export

The Penpot deployment and native editing workflow are documented in
`PENPOT_SETUP_AND_USAGE.md`. `fig1_penpot_build.js` records the programmatic origin
of the native layer hierarchy; final optical adjustments were made in Penpot.
Portable publication exports are normalized with:

```bash
uv run python figures/penpot/normalize_exports.py
```

The selected SVG is self-contained, the PDF is vector, and the inspection PNG is
3,600 × 2,240 pixels. The physical export size is 180 mm × 112 mm.
