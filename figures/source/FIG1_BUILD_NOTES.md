# Figure 1 build notes

> **Current master:** the final art-directed Figure 1 is the native Penpot file
> `figures/penpot/SolvAI_Figure1.penpot`. The SVGWrite workflow described below is
> retained as the deterministic scientific-asset generator and design predecessor,
> not as the final layout source. Penpot deployment, MCP editing, variants, exports
> and external-asset licensing are documented in
> `figures/penpot/PENPOT_SETUP_AND_USAGE.md`.

## Reproduction

From the repository root:

```bash
uv sync --extra dev
uv run python figures/source/fig1_build.py
```

The script renders three compositions at 180-mm publication width. Alternatives A
and B are written to `figures/alternatives/`; the selected balanced composition is
written as `figures/main/fig1_variant_C_balanced.*`, copied to
`figures/main/fig1_concept.{svg,pdf,png}`, and mirrored into `paper/figures/main/`.
Each variant also has a 300-dpi, 180-mm print-review PNG. `make figures` invokes the
same builder before creating the remaining paper figures.

## Scientific graphics toolchain

| Tool | Frozen version | Role |
|---|---:|---|
| RDKit | 2026.3.5 | Canonical molecular structures, 2D SVG depiction, ETKDGv3 conformers, MMFF optimisation and van der Waals radii |
| Indigo Toolkit | 1.46.0 | Independent 2D depiction comparison for the representative solute |
| Open Babel | 3.1.0 | Chemically valid 3D solute and water coordinate generation for the solvent-shell comparison |
| Packmol | 20.010 | Deterministic packing of the real solute--water configuration used in alternative B |
| PyMOL open source | 3.2.0a | Rendering comparison for the packed solute--water configuration |
| Matplotlib | 3.11.1 | Vector rendering of the measured lambda response, conformers, cavity contours and mathematical output label |
| pandas / pyarrow | 3.0.5 / 25.0.1 | Reading the frozen PIMD2 response table |
| SVGWrite | 1.4.3 | Exact-grid SVG composition and labelled response-coordinate artwork |
| CairoSVG | 2.8.2 | Vector PDF export and high-resolution PNG inspection render |
| Inkscape CLI | 1.4.2 | Independent SVG/PDF rendering and optical inspection |

All content is deterministic. Molecule/conformer generation uses fixed seeds, and
the exact source paths and generated hashes are recorded in
`fig1_assets/asset_manifest.json`.

The architecture is drawn directly in SVGWrite rather than with a generic neural
network glyph. Cairo's generated PDF timestamp is normalized to make repeat builds
byte-for-byte reproducible. Inkscape was used as an independent rendering check,
but CairoSVG remains the release exporter because it preserves the embedded source
SVGs as vectors whereas Inkscape 1.4.2 rasterized them during PDF import.

## Art-direction calibration

The redesign was calibrated against recent molecular-science and molecular-ML
figures from Nature Communications and Nature Machine Intelligence, including
DOIs `10.1038/s41467-025-57688-8`, `10.1038/s41467-024-55320-9`,
`10.1038/s41467-024-47999-7`, `10.1038/s42256-024-00856-0`,
`10.1038/s42256-025-01031-9` and `10.1038/s42256-024-00843-5`. The useful common
grammar was white space, direct molecular evidence, short panel-level claims,
thin connectors and low decorative overhead. No artwork was copied.

## Renderer and composition decisions

RDKit and Indigo depictions of N-methylacetamide are both preserved. Indigo's
explicit terminal methyl labels are clear in isolation, but RDKit's skeletal form
is more legible at the final panel size and was selected. Packmol generated a
72-water configuration from Open Babel coordinates (`water_shell.packmol.inp` and
`nma_water_shell.pdb`). Both an editable vector projection and a PyMOL render are
preserved. The molecular-science-forward composition uses the vector projection;
the PyMOL raster was rejected from the final figure because it introduced visual
weight and raster content without clarifying a retained SolvAI coordinate.

Three layouts were compared at 180 mm:

- **A, Minimal Nature-style:** strongest whitespace and fastest global reading,
  but it compresses the source chemistry too aggressively.
- **B, Molecular-science-forward:** foregrounds the real packed solvent shell and
  conformers, but overweights an illustrative configuration and the non-retained
  lambda diagnostic.
- **C, Balanced (selected):** retains four grounded source vignettes, replaces the
  tabular prior block with a four-family coordinate atlas, and gives deployment the
  largest uninterrupted visual path. It best preserves the training/deployment
  distinction without making the model look like either a software workflow or an
  explicit-solvent method.

## Scientific source assets

### Molecular structure

`n_methylacetamide_rdkit.svg` is an RDKit-authored depiction of
N-methylacetamide (`CNC(C)=O`), a neutral molecule in the ARROW reference set. The
cropping step changes only the SVG viewport; it does not redraw molecular bonds or
atom labels. `n_methylacetamide_rdkit_canvas.svg` preserves the original RDKit canvas
used for the cavity overlay.

### Calculated-solvation vignette

`continuum_cavity.svg` combines the RDKit molecular depiction with nested contours
computed from the molecule's RDKit atom coordinates and periodic-table van der Waals
radii. It is deliberately labelled a van der Waals cavity: it is a scientific
schematic of the molecular envelope used to explain continuum response, not a
calculated electrostatic potential or an SMD/COSMOtherm output.

### Empirical physicochemical coordinates

`abraham_axes.svg` shows the five named Abraham coordinates E, S, A, B and L. The
five spokes have equal length and encode categories only; they do not imply measured
magnitudes for the example molecule. The physical meanings are excess molar
refraction, dipolarity/polarizability, hydrogen-bond acidity, hydrogen-bond basicity
and hexadecane--air partition response.

### Conformational response

`dimethoxyethane_conformers.svg` shows three real conformers of benchmark
1,2-dimethoxyethane (`COCCOC`). Thirty conformers are generated with RDKit ETKDGv3
using seed `0x5A17`, optimised with MMFF, and three geometrically distinct low-energy
states are selected and aligned. Their coordinates are preserved in
`dimethoxyethane_selected_conformers.sdf`; conformer identifiers and relative MMFF
energies are recorded in `dimethoxyethane_conformer_metadata.json`. These geometries
illustrate conformational response; they are not substituted for the ConfSolv source
calculations.

### Explicit-solvent comparison asset

`nma_openbabel.pdb` and `water_openbabel.pdb` are the frozen Open Babel coordinate
inputs. Packmol 20.010 places 72 waters using the fixed seed in
`water_shell.packmol.inp`; the resulting coordinates are
`nma_water_shell.pdb`. `nma_water_shell_vector.svg` is a deterministic orthographic
projection of the 18 nearest packed waters, while `water_shell.pml` produces the
separate PyMOL comparison render. The exact regeneration commands are:

```bash
cd figures/source/fig1_assets
packmol < water_shell.packmol.inp
uv run pymol -cq water_shell.pml
```

The packed configuration is illustrative source artwork, not an equilibrated
trajectory or an input to SolvAI. It appears only in rejected composition B.

### Alchemical-response diagnostic

`nmethylacetamide_lambda_response.svg` plots the actual N-methylacetamide PIMD2
observations at lambda = 0.1, 0.5 and 0.9 from
`results/ablations/pimd2_multilambda_teacher.parquet`. The total, Coulomb and van der
Waals values are copied to `nmethylacetamide_lambda_response.csv`. The plot uses
straight segments between measured states and no invented smoothing. Its dashed
annotation states **explored, not retained**: no lambda, PIMD2, PIMD4, PIMD8 or NQE
coordinate occurs in the released 15-prior vector.

### The 15 response priors

`response_priors_15.svg` is generated from an explicit 15-entry specification
cross-checked against `paper/supplementary/tables/response_priors.tex` and
`data/manifests/final_training_config.json`. The entries are grouped as:

- two calculated-solvation responses: COSMOtherm water and SMD(water);
- five empirical polarity/hydrogen-bonding coordinates: Abraham E, S, A, B and L;
- two corrected water-model responses: OpenFF explicit water and GBn2 implicit
  solvent, each including its learned experimental residual;
- six ConfSolv conformational summaries: gas, solution and hydration conformer
  corrections; conformer solvation-energy spread; mean conformer solvent response;
  and conformer response spread.

They are named scalar targets, not hidden or latent embeddings.

## Architecture represented

Stage 1 contains the six released structure-to-response artifacts: two CHEMELEON
D-MPNNs, three ExtraTrees teachers and one LightGBM teacher. Stage 2 combines 2,048
Morgan bits and 217 RDKit descriptors (2,265 deterministic structure features) with
the 15 predicted responses. Experimental hydration labels supervise the three-member
ExtraTrees endpoint ensemble; each member contains 360 trees.

At deployment, the frozen response surrogates still run, but their original physical
calculations and source tables do not. The only user input is SMILES. No MD, PIMD or
probe calculation runs at inference. ARROW/PIMD8 is drawn on a separate magenta rule
as an accuracy comparator (approximately 0.205 kcal mol-1), with no arrow into the
training or inference graph.

## Visual system

- orange: source physical/physicochemical response information generated during
  teacher-data construction;
- blue: learned structure-to-response surrogates and their named outputs;
- teal: experimental endpoint supervision and deployable SolvAI inference;
- dark neutral: structure and deterministic molecular information;
- magenta: ARROW/PIMD8 accuracy reference only.

The master is a self-contained SVG with embedded vector source assets. It uses an
180-mm by 105-mm canvas, no gradients, shadows, clip art or raster scientific
content. The PNG is a 3,600-pixel-wide inspection derivative; the PDF remains vector.
