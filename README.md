# SolvAI

**SolvAI learns reusable solvent-response coordinates from calculated, empirical and
corrected solvation data and predicts hydration free energy directly from molecular
structure—without running simulation at inference.**

![SolvAI concept](paper/figures/main/fig1_concept.svg)

The released system maps one SMILES string to one hydration free energy. Its response
surrogates were trained on benchmark-disjoint quantum-continuum, alchemical,
empirical and conformational data; the expensive source calculations are not rerun
for a query. PIMD-derived features were tested but are not present in the final model.

## Confirmatory result

On the 85-solute neutral-hydration reference set introduced with ARROW:

| Method | Simulation at inference? | MAE (kcal/mol) |
|---|---:|---:|
| Classical ARROW | yes | 0.785 |
| ARROW/PIMD8 | yes | 0.205 |
| Matched structure-only endpoint | no | 0.303 |
| SolvAI, fixed five-fold OOF | **no** | **0.202** |
| SolvAI, five complete partitions | **no** | **0.207 ± 0.004** |
| SolvAI, no ARROW labels in training | **no** | **0.257** |

The matched endpoint uses exactly the same experimental labels, descriptors,
ExtraTrees architecture, weights, folds and seeds; only the 15 response priors are
removed. The paired OOF improvement is −0.101 kcal/mol (95% bootstrap interval,
−0.215 to −0.020). Shuffled priors do not improve the endpoint, and the advantage
survives global family, scaffold, molecular-cluster and nearest-neighbour exclusions.

In a prospectively frozen external molecule-disjoint cohort, the same matched
comparison lowers MAE from 1.532 to 1.153 kcal/mol (N=220). The advantage also
survives in 97 molecules absent from all six response-teacher source tables (2.138 to
1.536 kcal/mol). These broader molecules establish transfer of the response-layer
advantage, not PIMD8-level absolute accuracy outside ARROW-85.

The supported conclusion is PIMD8-level accuracy on this reference chemistry, not a
general sub-0.20 claim. Global family and scaffold separation
remain harder at 0.468 and 0.376 kcal/mol, respectively.

## Install and predict

Python 3.11, [uv](https://docs.astral.sh/uv/) and Git LFS are required.

```bash
git lfs install
git clone https://github.com/Scientific-Computing-Lab/SolvAI.git
cd SolvAI
make setup
uv run solvai predict 'CCO'
```

The command returns the ensemble-mean hydration free energy and ensemble spread in
kcal/mol:

```text
CCO    -5.012566    0.004714
```

The API is equally small:

```python
from solv_ai import predict_smiles

prediction, spread = predict_smiles(["CCO", "c1ccccc1"])
```

On the release CPU host, the packaged artifact has a 15.29 s warm single-molecule
median and processes a batch of 32 in 15.82 s (0.494 s per molecule). Startup of the
two response D-MPNNs dominates single-query latency; no simulation is performed.

## Reproduce the paper

```bash
make test && make verify && make figures && make paper
```

This quick path uses frozen, hash-verified artifacts to recompute predictions,
metrics, tables, figures and PDFs. It does not rerun physical calculations or model
training. The preregistered confirmation protocol is in
[`release/CONFIRMATORY_FREEZE.md`](release/CONFIRMATORY_FREEZE.md), with results in
[`reports/CONFIRMATORY_ANALYSIS.md`](reports/CONFIRMATORY_ANALYSIS.md).
The prospectively frozen external protocol and report are
[`release/TIER_A_EXTERNAL_VALIDATION_FREEZE.md`](release/TIER_A_EXTERNAL_VALIDATION_FREEZE.md)
and [`reports/TIER_A_EXTERNAL_VALIDATION.md`](reports/TIER_A_EXTERNAL_VALIDATION.md).

The compiled [manuscript](paper/main.pdf),
[Supplementary Information](paper/supplementary/supplementary.pdf), including five
supporting figures, and machine-readable Supplementary Data are included. See
[`repro/QUICK_REPRODUCTION.md`](repro/QUICK_REPRODUCTION.md),
[`repro/FULL_REPRODUCTION.md`](repro/FULL_REPRODUCTION.md) and
[`repro/DATA_PROVENANCE.md`](repro/DATA_PROVENANCE.md).

## Scientific safeguards

- Exact and standardized benchmark equivalents are absent from all supervised
  external training sources used by the confirmatory model.
- Every reported accuracy value is held out; the all-data deployment refit is never
  used as evidence.
- Shuffled-prior, global chemical-separation and zero-ARROW-label controls are
  included molecule by molecule.
- Tier-A eligibility was frozen before evaluation; all 220 rows are endpoint-disjoint
  and the strict 97-molecule subset is also response-source-disjoint.
- Inference requires no experimental target, family/scaffold label, MD, PIMD, ARROW
  trajectory, probe or routing policy.
- The released artifact contains no retained PIMD-trained feature.
- The returned ensemble spread is not a calibrated applicability or reliability
  score; the validated scope is neutral small-molecule hydration chemistry.

## Repository map

- `solv_ai/` — SMILES-only inference and metric code
- `models/final/` — standardized-exclusion response surrogates and endpoint ensemble
- `results/confirmatory/` — preregistered predictions, comparisons and statistics
- `audits/confirmatory/` — identity, similarity and refit audits
- `paper/` — Nature Communications manuscript and Supplementary files
- `repro/` — quick/full reproduction and data provenance

Citation metadata are provided in `CITATION.cff`. Code is MIT licensed; external
datasets retain the terms listed in the provenance record.
