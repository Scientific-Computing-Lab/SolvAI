# Round 1 — Reviewer D: statistics and skeptical generalist

## Critical issues

None after the repeated-split result was placed in the headline figure. The
paper no longer treats 0.20 kcal mol−1 as a statistical decision boundary.

## Important issues

- State the 0.197 result as a strict fixed-partition point estimate and report
  the five-partition centre, 0.204 ± 0.005, in the same visual field.
- The difference from PIMD8 is not resolved by paired resampling. Use
  “PIMD8-level” rather than “better than PIMD8”.
- Forty-two of 85 individual absolute errors decrease although the aggregate
  MAE improves. Explain that the gain is magnitude-weighted rather than imply a
  majority vote.
- Small chemical families must display sample size and should not be precisely
  ranked.

## Optional polish

- Keep bootstrap intervals in Extended Data; the split-repeat points are more
  immediately informative in the main figure.

## Revision made

Figure 2 shows all repeat values without connecting them. Figure 4 and Extended
Data Figure 6 show molecule-level family errors and n. The paired-bootstrap
intervals and fixed-versus-nested intervals are consolidated in Extended Data
Figure 7.
