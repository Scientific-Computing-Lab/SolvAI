# SolvAI model card

SolvAI predicts neutral-molecule hydration free energy from a SMILES string. The
released stack computes deterministic molecular descriptors and 15 response priors
with structure surrogates, then applies an ensemble of three ExtraTrees endpoint
models. No MD, PIMD, ARROW trajectory, probe calculation or experimental lookup is
performed at inference.

RDKit parses each query and emits canonical isomeric SMILES before feature
generation. This is not a claim of invariance across tautomers, protonation states or
salt forms.

## Validated performance

- Confirmatory fixed five-fold OOF MAE: 0.20223 kcal/mol
- Matched no-prior endpoint MAE: 0.30335 kcal/mol
- Five-partition mean: 0.20737 ± 0.00444 kcal/mol
- Zero-ARROW-label transfer MAE: 0.25694 kcal/mol
- Tier-A endpoint-disjoint MAE (N=220): 1.15255 kcal/mol
- Tier-A strict source-disjoint MAE (N=97): 1.53560 kcal/mol
- Global family / scaffold MAE: 0.46779 / 0.37619 kcal/mol

The reference domain is neutral organic hydration chemistry. The model is not
validated for ions, salts, metals, proteins or broad chemical extrapolation. The
fixed point estimate is on the ARROW/PIMD8 accuracy scale, but robust sub-0.20
performance and superiority over PIMD8 are not claimed. The returned ensemble spread
is not a calibrated per-query applicability or reliability score.

## Training and inference boundary

Physical calculations supplied training targets for the response surrogates. The
deployed model receives structure only. PIMD-derived candidate features were tested
and not retained; PIMD8 is used solely as an accuracy comparator.
