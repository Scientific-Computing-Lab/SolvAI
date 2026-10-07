# SolvAI journal manuscript

`main.tex` is the single compilation root. Its output, `main.pdf`, always contains
the main text, references and all Supplementary Information in that order, with
continuous page numbering and no reviewer line numbers. The supplement shares the main document's
font, margins and paragraph settings. Its S-numbered methods, notes, tables and
figures retain their existing identifiers.

From the repository root, with LaTeX and latexmk installed:

```sh
make journal-paper
```

Alternatively, from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Edit the main text in `main.tex` and the supplement in
`supplementary/supplementary.tex`. The latter is an input fragment, not a separate
document. Do not concatenate separately compiled PDFs: a single LaTeX build keeps
pagination, bookmarks and cross-references consistent. There is no independent
supplement PDF to maintain.

The publication build includes eight main figures and twelve supplementary figures.
Quantitative plots are vector PDF/SVG. Figures 1--2 and Supplementary
Figures 8--10 combine conceptual raster motifs with editable vector labels and
connectors; they are not quantitative data or molecular structures. Components and
composition scripts are retained. These sources do not contain author-facing submission
notes, editorial reviews or internal design prompts. Machine-readable supplementary
data remain in `supplementary_data/` and are not replaced by the PDF tables.
The extended endpoint study is in `endpoint_results.tex`, `endpoint_methods.tex`
and `supplementary/endpoint_methods.tex`; its source outputs are in
`../results/endpoint_models_20261006/`. See `../repro/endpoint_models/README.md`
for the distinction between reproducible frozen results and licensed training inputs.
