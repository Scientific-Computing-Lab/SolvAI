# Quick reproduction

This path uses frozen, hash-verified artifacts. It does not rerun quantum
chemistry, molecular dynamics, PIMD or model training.

Requirements: Linux, Python 3.11, Git LFS, `uv`, a LaTeX installation with
`latexmk`, and `qpdf`; allow approximately 3 GB for the Python environment and
400 MB for the repository and artifacts. A GPU is not required.

```bash
git clone https://github.com/Scientific-Computing-Lab/SolvAI.git
cd SolvAI
make setup
make test
make verify
uv run solvai predict 'CCO'
make figures
make paper
```

`make verify` recomputes the paper metrics from molecule-level predictions,
rechecks the available benchmark-disjoint source tables and model schema, and
verifies artifact hashes and frozen SMILES-only predictions. `make figures`
regenerates every main and Extended Data panel. `make paper` compiles the main
manuscript and Supplementary Information.

Expected ethanol output (prediction and ensemble spread, kcal/mol):

```text
CCO    -5.020248    0.005697
```

Small floating-point differences are accepted only within the
$2\times10^{-6}$ kcal/mol tolerance recorded in the artifact manifest.
