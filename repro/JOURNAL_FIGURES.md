# Journal figures

The current main manuscript uses eight figures in `paper/figures/journal/`.
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

## Scientific provenance

All plotted numbers come from the frozen confirmatory and journal result tables.
`paper/figures/journal/sources.json` records input and builder SHA-256 hashes.
The molecular illustration uses retained N-methylacetamide connectivity and
coordinates. Conformers use the retained 1,2-dimethoxyethane SDF. Neither solvent
shell placement nor schematic energy levels is represented as simulation data.
Quantitative plots contain no embedded rasters. Figures 1--2 and Supplementary
Figures 8--10 combine selected Azure-generated conceptual motifs
with native SVG labels, arrows and scientific annotations. These are explicitly
non-evidential illustrations, not molecular structures or measured data. Selected
components are in `paper/figures/journal/teaser_components/` and
`paper/supplementary/figures/endpoint_components/`; the composition
script rebuilds the exports without an image-generation API or credentials.

The paper-wide plotting style is in `scripts/journal_style.py`; opening-figure
composition is in `scripts/journal_teaser.py`. Native SVG labels, molecule counts,
operators and arrows carry the scientific meaning; translucent objects depict
information processing rather than physical fields or exact network layer counts.

The main sequence is: overview; learning and evaluation boundaries; matched
evidence; molecular-model comparisons; chemical generalization; source-composition
reversal; cross-endpoint response ablation; labelled size transfer and unlabelled
peptide series. Captions in the manuscript define every panel and interval. The
supplement retains all original exploratory and confirmatory diagnostics and adds
the graph-model selection/refit protocol, conformer-target interpretation, nine
endpoint mechanisms, native interval coverage, and duration/capacity challenges.
`results/endpoint_models_20261006/` contains the corresponding numeric sources and
the corrected fixed-tree MoLFormer comparison. Original archived outputs remain
unchanged and are distinguished from corrected scientific comparisons.

Figures preserve all nine model conditions, all four similarity thresholds/bins,
all five repeat partitions, all source-block contrasts and the external strict
subset. Strict-97 is nested in external-220. The computed-only control retains ten
descriptors by removing five Abraham values and two residual corrections, not two
physical-energy descriptors.
