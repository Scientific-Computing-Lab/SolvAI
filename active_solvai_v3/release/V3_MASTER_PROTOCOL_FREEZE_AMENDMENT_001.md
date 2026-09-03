# V3 master protocol freeze amendment 001

**Committed before any Gate-1 model fitting or held-out scoring:** 2026-09-03  
**Information known:** inherited aggregate results, raw trajectory inventory and
the frozen twelve-molecule feature/time-series tables; no v3 target-model result.

## Clarification

The master freeze says that structure- and SolvAI-conditioned models add their
fixed molecule-level coordinates to lambda/protocol features. An additive
linear model would allow those coordinates to shift average difficulty between
molecules but could not represent a molecule-conditioned *ranking of lambda
windows*. That would fail to implement the registered scientific question.

Accordingly, the predeclared linear feature map is clarified as follows:

- every molecule-level structure or SolvAI coordinate enters as a main effect;
- it also enters multiplied by `(lambda - 0.5)` and `(lambda - 0.5)^2`;
- no outcome-dependent interaction selection is performed;
- all interactions use the same training-only imputation, standardization,
  grouped validation and ridge-alpha grid already frozen;
- the lambda/protocol and generic observed-only models are unchanged.

The molecule-shuffled control applies the same interaction map after jointly
permuting the complete response vector among outer-training molecules. The
held-out molecule's response vector is not inserted into any training row.

This clarification changes neither the target, data, folds, penalties, metrics
nor pass rules. It prevents an implementation artifact from making the
molecule-conditioned hypothesis algebraically untestable.
