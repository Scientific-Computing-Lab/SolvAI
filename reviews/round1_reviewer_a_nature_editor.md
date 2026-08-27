# Round 1 — Reviewer A: Nature editor

## Critical issues

1. The broad advance must precede model details. The initial framing risked
   reading as a small benchmark improvement. Lead with the change in how
   simulation is used: reusable training supervision instead of per-query
   computation.
2. The title and abstract must not imply that the selected artifact directly
   uses PIMD labels. It does not; it reaches the PIMD8 reference accuracy using
   other solvent-response teachers.

## Important issues

- Keep the 85-solute scope visible without calling it a community-standard
  benchmark.
- Treat the 0.197 result and 0.204 repeat mean together.
- Reduce the main narrative to aligned response, headline result, transfer and
  frontier; move campaign history to Supplementary Information.

## Optional polish

- Make Figure 1 readable before its caption.
- End with the reusable-supervision concept, not a performance slogan.

## Revision made

The abstract and opening now lead with amortized physical supervision. The text
states that zero PIMD8 labels are present in the selected artifact, reports the
five-repeat mean beside the headline point estimate, and limits the principal
claim to the ARROW reference chemistry.
