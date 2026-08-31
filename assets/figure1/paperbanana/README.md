# PaperBanana design record

These files record the PaperBanana portion of the Figure 1 workflow:

- `brief.md` and `caption.txt` are the source-grounded inputs.
- `planner_result.json` contains the PaperBanana Planner's full A–F design
  specification.
- `critic_round1.json` reviews the deterministic vector draft.
- `critic_round2.json` reviews the first hybrid SVG render with embedded ImageGen
  assets.
- `critic_round3.json` reviews the asset-rich hybrid render and supplies the
  final scientific-consistency pass.
- `panel_c_revision_brief.md` records the author-supplied scientific correction
  that supersedes the earlier varied-molecule Panel C concept.
- `critic_round4.json` verifies the revised one-molecule, six-surrogate Panel C
  and its exact `1 | 5 | 1 | 1 | 1 | 6` concatenation.
- `critic_round5.json` verifies that panels D/E reuse the corrected wedge
  vocabulary and that both E-stage models expose one clean inference port.
- `critic_round6.json` verifies the flat light v2 wedges and the explicit D/E
  supervision and merge-arrow endpoints.
- `critic_round7.json` verifies the richer learned-forest endpoint asset and the
  white-only sampling enclosure in panel A.
- `critic_round8.json` verifies the corrected endpoint topology: all three input
  blocks concatenate before one 2,280-D representation is evaluated by three
  full-feature ExtraTrees ensembles and their predictions are averaged.

The vector compositor lives in `scripts/make_figure1_overview.py`. This adapter
preserves PaperBanana's Planner→Visualizer→Critic idea while replacing its
whole-canvas raster visualizer with a deterministic SVG compositor. Twenty-one
curated, text-free visual assets are embedded as PNG image elements; layout,
labels, connectors and the quantitative evidence panel stay editable and
source-grounded. The Critic's scientific corrections are applied selectively;
stylistic suggestions that conflict with the requested asset-rich design are
not treated as requirements.
