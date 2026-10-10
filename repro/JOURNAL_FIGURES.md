# Journal figures

The current main manuscript uses seven figures in `paper/figures/journal/`.
Every figure has an editable SVG and a publication PDF. The Supplementary Information
uses twelve figures in `paper/supplementary/figures/` and is automatically included in
the same `paper/main.pdf` by `make journal-paper`.

## Build

Run from the repository root in a Python environment with NumPy, Pandas, PyArrow,
Matplotlib, RDKit and CairoSVG. The Cairo library must be available to CairoSVG.

```sh
python scripts/make_figures.py
```

For the journal display sequence from retained results and illustration components:

```sh
python scripts/make_endpoint_figures.py
python scripts/make_local_journal_figures.py
python scripts/make_endpoint_concepts.py
```

The optional `--preview-dir PATH` creates lightweight raster previews **outside**
the manuscript's publication dependencies. On an Apple Silicon Mac with Homebrew
Cairo, set `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` if needed. Arial is used
for the local publication build. The scripts require no GPU or remote computation.

For an art-only refresh of the quantitative displays, without rewriting result
tables or the conceptual figures, run:

```sh
python scripts/make_endpoint_figures.py --plots-only
python scripts/make_local_journal_figures.py --only matched modern transfer composition
python scripts/make_figures.py --quantitative-supplement-only
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/main.tex
```

This path covers main Figures 3–8 and Supplementary Figures 1, 3–5, 11–12.
`scripts/journal_style.py` sets their common scientific plotting language:
muted model and response colors, thin axes and intervals, light reference grids,
white backgrounds, consistent units and panel typography. Exact points, intervals,
sample sizes and evaluation partitions remain defined by the retained numeric
sources; decorative shading and redundant per-point value labels have been
removed where they obscured comparisons. Illustrative conceptual panels remain
separate from data-derived plot panels.

## Scientific provenance

All plotted numbers come from the frozen confirmatory and journal result tables.
`paper/figures/journal/sources.json` records input and builder SHA-256 hashes.
The molecular illustration uses retained N-methylacetamide connectivity and
coordinates. Conformers use the retained 1,2-dimethoxyethane SDF. Neither solvent
shell placement nor schematic energy levels is represented as simulation data.
Quantitative plots and all mechanism diagrams are native vectors. The merged
Figure 1 uses six Azure-generated conceptual source icons, a solvent halo, an
architecture-neutral endpoint lens and six paper-texture samples. The textures
are clipped inside native vector ribbons; their SVG embedding uses one image
per ribbon to avoid tiled seams. Ribbon geometry, verified molecular structures,
labels, connections and coordinate counts are authored, not generated. The halo
and lens have exact-white backgrounds after deterministic processing with
`scripts/prepare_teaser_components.py`; the raw Azure output and both sanitized
provenance records remain in `paper/figures/journal/teaser_components/`.
Figure 1c redraws the exact 15 coordinates as native colored tiles. Its NMA,
fingerprint-bit and RDKit-feature glyphs illustrate types of structural input;
the bit pattern and bar heights are not a measured feature vector. Figure 1d
contains exactly 85 ARROW tiles and 220 external tiles, the latter partitioned
into 123 source-exposed and 97 source-disjoint tiles. Tile placement is a schematic
record layout, not a molecular ordering.
Individual vector mechanism panels are in
`paper/supplementary/figures/endpoint_components/`. Rebuilding the exports needs
neither an image-generation API nor credentials.

The paper-wide plotting style is in `scripts/journal_style.py`; opening-figure
composition is in `scripts/journal_teaser.py`. `scripts/journal_diagram.py` defines
native function ports, tree geometry, concatenation and addition. Diagram connectors
join the intended glyph or model boundaries. Colored streams depict information reuse, not
physical fields; generic hydration functions do not imply a particular architecture.

The main sequence is: merged physical/source/endpoint/evaluation overview; matched
molecule-alignment evidence; cross-endpoint response comparison; molecular-model
comparisons; chemical generalization; source-composition reversal; labelled size
transfer and unlabelled peptide series. Captions in the manuscript define every panel
and interval. The supplement retains exploratory and confirmatory diagnostics and adds
the graph-model selection/refit protocol, conformer-target interpretation, nine
endpoint mechanisms, native interval coverage, and duration/capacity challenges.
Its chemical-family error figure now omits the global-separation plot already in
the main chemical-generalization figure.
`results/endpoint_models_20261006/` contains the corresponding numeric sources and
the corrected fixed-tree MoLFormer comparison. Original archived outputs remain
unchanged and are distinguished from corrected scientific comparisons.

Figures preserve all nine model conditions, all four similarity thresholds/bins,
all five repeat partitions, all source-block contrasts and the external strict
subset. Strict-97 is nested in external-220. The computed-only control retains ten
descriptors by removing five Abraham values and two residual corrections, not two
physical-energy descriptors.
