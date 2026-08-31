# Figure 1 generated asset register

Figure 1 is a hybrid scientific illustration. The layout, text, arrows and
quantitative panel are created deterministically by
`scripts/make_figure1_overview.py`. Twenty-one text-free raster assets provide the
molecules and pipeline components and are embedded directly into the final SVG
as base64 PNG data.
The current composition contains 30 embedded image instances drawn from 21
selected unique assets.

PaperBanana is used as the design controller: its Planner produced the A–F
composition specification and its Critic reviews rendered SVG iterations against
the scientific brief and caption. The standard PaperBanana image visualizer is
not used for the whole figure because it returns a flattened raster image.

## Selected assets

- `acetamide_3d_v2.png`: isolated acetamide, exact requested connectivity
  `CH3–C(=O)–NH2`, on a uniform white background.
- `acetamide_hydrated_v1.png`: composite of the selected acetamide and a ring of
  discrete water molecules; used only to illustrate the physical hydration
  calculation in panel A.
- `response_manifold_v1.png`: six text-free amber response surfaces converging
  into an ordered teal channel bundle; used as a low-opacity conceptual layer in
  panel B, not as quantitative data.
- `sampling_frames_v1.png`: three per-solute explicit-water sampling snapshots.
- `cosmo_continuum_v1.png`, `abraham_axes_v1.png`, `openff_explicit_v1.png`,
  `gbn2_implicit_v1.png`, `smd_continuum_v1.png` and
  `confsolv_conformers_v1.png`: the six response-teacher vignettes.
- `panelc_query_acetophenone_v1.png`: the single neutral-organic query used once
  at the entrance to all six Panel C branches.
- `model_icon_dmpnn_azure_v1.png`, `model_icon_extratrees_azure_v1.png` and
  `model_icon_lightgbm_azure_v1.png`: high-contrast Azure-generated model-class
  icons used for all six Panel C branches. Reuse makes the repeated architecture
  legible and avoids the former flattened miniature scenes.
- `response_channels_15_v1.png`: reusable response-vector asset in panels D/E;
  Panel C's exact family-owned response cells are drawn as native SVG geometry.
- `fingerprint_v1.png` and `descriptors_v1.png`: endpoint structure-feature
  assets.
- `endpoint_extratrees_triplet_azure_v1.png`: one shared input that branches to
  three parallel ExtraTrees learners and reconverges at one mean-prediction
  node; used for the trained endpoint ensemble in panels D and E. The topology
  prevents the three input feature blocks from being misread as three separate
  model inputs.
- `frozen_surrogate_wedge_flat_v2.png`: a pale, flat six-lane surrogate-bank
  wedge with exactly one inlet and one outlet; used in panel E.
- `lock_icon_flat_azure_v1.png`: a clean, flat padlock generated with the Azure
  `gpt-image-1.5` deployment. The compositor removes its white background and
  recolours the same mask blue or green for the two frozen-model locations.
- `experimental_dg_strip_v1.png`: text-free experimental endpoint supervision.

## Supporting/draft assets

- `hydration_shell_v1.png`: the water-shell source used to compose the hydrated
  acetamide asset.
- `acetamide_3d_v1.png`: rejected checkerboard-background draft; retained only
  for provenance and not embedded in the final figure.
- `molecular_graph_v1.png`, `molecular_graph_phenol_v1.png`,
  `molecular_graph_acetone_v1.png` and `surrogate_module_v1.png`: superseded
  Panel C assets retained for provenance but no longer embedded.
- `panelc_combisolv_dmpnn_v1.png`, `panelc_smd_dmpnn_v1.png`,
  `panelc_abraham_extratrees_v1.png`, `panelc_openff_extratrees_v1.png`,
  `panelc_gbn2_extratrees_v1.png` and `panelc_confsolv_lightgbm_v1.png`:
  superseded low-contrast Panel C mechanism scenes retained for provenance.
- `endpoint_ensemble_v1.png` and `frozen_model_v1.png`: superseded network/cable
  motifs retained for provenance but no longer embedded in panels D/E.
- `endpoint_extratrees_wedge_v1.png` and `frozen_surrogate_wedge_v1.png`:
  superseded dark three-dimensional wedge drafts retained for provenance.
- `endpoint_extratrees_wedge_flat_v2.png`: superseded single-tree-per-wedge draft;
  retained for provenance after the learned-forest v3 replacement.
- `endpoint_extratrees_forest_flat_v3.png`: superseded three-row endpoint forest;
  retained for provenance because its rows could be misread as corresponding to
  the three feature blocks.

All scientific labels and numerical values remain outside the generated assets.
The final SVG must contain exactly 30 embedded raster image elements.
Its editable SVG text uses the same Latin Modern Roman and Latin Modern Mono
families as the LaTeX manuscript; the regular, bold and monospaced OpenType files
are embedded in the SVG so the typography does not depend on fonts installed on
the viewing computer. The generator also rejects text overlaps and text that
extends outside a panel before writing any output.

## ImageGen prompt contract

The assets were generated with Codex's built-in ImageGen mode, one component at
a time. Every prompt used the same production contract: isolated text-free
scientific component; clean Nature-style editorial rendering; uniform pure-white
background; no border, caption, panel letter, watermark, checkerboard or drop
shadow; generous whitespace; chemically and geometrically legible at small size.
The original prompt subjects were the exact concepts listed above: acetamide in gas
and hydrated contexts, six distinct response teachers, molecular examples, a frozen
structure-to-response surrogate, response channels, Morgan fingerprint, RDKit
descriptor field, endpoint ensemble and experimental hydration-free-energy strip.

The Panel C revision used seven additional component prompts under the same contract:

- one accurate, isolated neutral acetophenone ball-and-stick query molecule;
- a text-free CHEMELEON/D-MPNN message-passing vignette for CombiSolv-QM;
- a distinct text-free CHEMELEON/D-MPNN message-passing vignette for SMD(water);
- three distinct text-free ExtraTrees forest/ensemble vignettes for Abraham,
  OpenFF and GBn2;
- one text-free sequential boosted-tree vignette for ConfSolv/LightGBM.

The D/E cleanup used two further prompts:

- a three-member ExtraTrees ensemble rendered as clean faceted decision-tree
  wedges converging to exactly one output hub;
- a single six-lane frozen-surrogate wedge with one flush inlet and one flush
  outlet and no exposed cables.

The final v2 prompts additionally required flat 2D vector-like rendering, white or
near-white interiors, thin outlines, no perspective, no extrusion, no bevels, no
shadows and no dark filled bodies. Merge nodes, labels and directional connectors
remain native SVG geometry.

The final Azure endpoint prompt required one common input node, three parallel
ExtraTrees forest learners, and one common mean-prediction node. It also required
varied-depth decision trees, crisp flat line art, a pure-white field, and no text,
boxes, shadows or dimensional styling. The compositor supplies the exact
`360 trees each`, seeds `11/29/47`, concatenation and averaging labels as SVG.

The final lock asset was generated separately with the Azure `gpt-image-1.5`
deployment as a centered, symmetric, flat padlock with a thin outline and one
keyhole on a pure-white field. The compositor converts it to an antialiased alpha
mask, so the published blue and green instances share identical geometry and no
generated background or shading remains visible.

The three Panel C model-class icons were also generated with Azure
`gpt-image-1.5`: a cyclic molecular message-passing graph for D-MPNN, three
decision trees converging to one output for ExtraTrees, and a left-to-right
additive tree sequence for LightGBM. Each prompt required high-contrast flat line
art, no text and no dimensional effects so the icons remain identifiable at final
print size.

Only the mechanism region of each surrogate image is cropped into the SVG, so no
generated molecule can contradict the single-query-molecule requirement. The exact
family names, surrogate types, branch ownership and `1 | 5 | 1 | 1 | 1 | 6 = 15`
accounting are deterministic SVG text and geometry, never generated pixels.
