# SolvAI model card

SolvAI predicts neutral-molecule hydration free energy from a SMILES string. The
released stack computes deterministic molecular descriptors and 15 response priors
with structure surrogates, then applies an ensemble of three ExtraTrees endpoint
models. No MD, PIMD, ARROW trajectory, probe calculation or experimental lookup is
performed at inference.

## Validated performance

- Confirmatory fixed five-fold OOF MAE: 0.20223 kcal/mol
- Matched no-prior endpoint MAE: 0.30335 kcal/mol
- Five-partition mean: 0.20737 ± 0.00444 kcal/mol
- Zero-ARROW-label transfer MAE: 0.25694 kcal/mol
- Global family / scaffold MAE: 0.46779 / 0.37619 kcal/mol

The reference domain is neutral organic hydration chemistry. The model is not
validated for ions, salts, metals, proteins or broad chemical extrapolation. The
fixed point estimate is on the ARROW/PIMD8 accuracy scale, but robust sub-0.20
performance and superiority over PIMD8 are not claimed.

## Training and inference boundary

Physical calculations supplied training targets for the response surrogates. The
deployed model receives structure only. PIMD-derived candidate features were tested
and not retained; PIMD8 is used solely as an accuracy comparator.
