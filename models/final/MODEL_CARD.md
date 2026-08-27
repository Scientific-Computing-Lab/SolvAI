# SolvAI model card

## Model summary

SolvAI is a structure-only model stack for predicting neutral-molecule hydration
free energy in water. Its public interface accepts SMILES and returns an
ensemble-mean prediction and ensemble spread in kcal/mol.

- Version: 1.0.0
- Input: one or more valid molecular SMILES strings
- Output: predicted hydration free energy, kcal/mol
- Inference: CPU-compatible; no molecular simulation
- Primary reference domain: neutral organic molecules represented by the
  ARROW 85-solute set and benchmark-disjoint teacher sources

## Architecture

The input is transformed into 2,265 deterministic RDKit/Morgan features. Six
bundled teacher artifacts predict 15 physical-response priors spanning quantum
continuum water response, Abraham solute parameters, alchemical and implicit
response, SMD(water), and ConfSolv conformer response. Three 360-tree
ExtraTrees regressors consume the resulting 2,280 features; their mean is the
reported prediction and their standard deviation is returned as model spread.

This is a model stack, not one monolithic neural network. The user-facing
operation remains one SMILES input to one hydration-free-energy output.

## Evaluation

| Evaluation | MAE (kcal/mol) |
|---|---:|
| Strict fixed five-fold OOF | 0.19705 |
| Nested feature selection OOF | 0.19931 |
| Five independent splits, fixed | 0.20375 ± 0.00493 |
| Five independent splits, nested | 0.20860 ± 0.00275 |
| Family holdout | 0.23957 |
| Scaffold holdout | 0.24128 |

The exact sub-0.20 point estimate is split-sensitive. The model should be
described as reaching PIMD8-level accuracy on this reference set, not as
robustly sub-0.20 or universally state of the art.

## Required and prohibited inputs

Required:

- parseable SMILES

Not required and not accepted:

- experimental hydration free energy;
- MD or PIMD trajectories;
- ARROW calculations or parameters;
- probe simulations;
- functional-group or scaffold labels;
- benchmark fold metadata.

All physical quantities used by the endpoint at deployment are predictions
from bundled structure-to-response models.

## Intended use

The model is intended for research predictions and method development on
neutral organic hydration chemistry. It can rank or estimate molecules whose
chemistry resembles the combined training domains. Ensemble spread is a useful
diagnostic but has not been calibrated as a guaranteed prediction interval.

## Limitations

- The primary evaluation has only 85 molecules.
- Family and scaffold extrapolation errors are about 0.24 kcal/mol.
- Amides, aromatics, ethers, acids and alkanes are the largest-error families
  on the frozen random OOF partition.
- Charged species, salts, mixtures, alternative protonation states and
  out-of-domain elements are not validated.
- A valid numerical output is not proof that a molecule lies in domain.
- The final selected artifact uses zero PIMD8 labels; sparse PIMD supervision
  was evaluated but did not improve held-out prediction.

## Reproducibility

`manifest.json` fixes the byte length and SHA-256 of every executable artifact
and records two end-to-end smoke predictions. Run:

```bash
make verify
```

The validation predictions come from models that did not receive each held-out
experimental target. The packaged deployment refit uses all 85 reference
labels only after evaluation and must not be used to recompute OOF accuracy.

## Ethics and responsibility

SolvAI is a research model. It does not replace experimental validation or
chemical safety assessment and should not be used as the sole basis for
clinical, environmental or safety-critical decisions.
