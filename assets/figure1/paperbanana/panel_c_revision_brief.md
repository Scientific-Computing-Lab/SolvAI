# Authoritative Panel C revision brief

This author-supplied correction supersedes the earlier PaperBanana suggestion to use
varied molecular structures in Panel C.

## Required scientific topology

One neutral organic query molecule feeds six independent frozen surrogate families:

1. CombiSolv-QM — CHEMELEON/D-MPNN — 1 coordinate.
2. Abraham — ExtraTrees — 5 coordinates.
3. OpenFF corrected response — ExtraTrees — 1 coordinate.
4. GBn2 corrected response — ExtraTrees — 1 coordinate.
5. SMD(water) — CHEMELEON/D-MPNN — 1 coordinate.
6. ConfSolv — LightGBM — 6 coordinates.

Each branch receives the same molecule and produces only its own coordinate group.
The six groups are concatenated exactly as `1 | 5 | 1 | 1 | 1 | 6 = 15` into one
molecule-aligned response vector.

## Required inference wording

“predicted per molecule by frozen structure→response surrogates · source calculations
are not run at inference”

The frozen surrogates are evaluated for every query molecule. The expensive source
calculations and simulations are what do not run at inference.

## Visualizer constraint

Image generation supplies only seven text-free raster subassets: one neutral query
molecule and six surrogate-mechanism vignettes. All family names, surrogate types,
counts, connectors, the concatenation rail and the 15-cell output remain native SVG.
The canvas background is pure white.
