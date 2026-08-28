# Skeptical novelty editorial check

**Perspective:** a *Nature Communications* editor familiar with molecular ML,
solvation modelling and physics-informed representations.

## Decision question

After the closest literature is acknowledged, is it clear why this study needed to
exist?

**Yes, with a bounded claim.** Prior work established direct hydration learning,
computed-to-experimental transfer, solvent-conditioned representations evaluated with
per-query physical calculations, and two-stage models that predict physical
descriptors from structure. Those precedents remove any defensible architectural
priority claim. They do not answer whether a response layer assembled across distinct
solvation formalisms contributes molecule-aligned information beyond an otherwise
identical experimentally supervised endpoint after all source calculations have been
removed from inference.

SolvAI answers that question with evidence that is unusual for a descriptor paper:

- a matched no-prior endpoint in which only the response layer changes;
- a shuffled-prior control that destroys molecular alignment without changing feature
  count;
- five complete partitions;
- family, scaffold, cluster and similarity separation applied to the complete
  endpoint-labelled pool;
- transfer with no ARROW labels; and
- negative high-fidelity and latent-response results that identify a mechanistic
  boundary.

The paper therefore supports a scientific claim about what physical information can
survive amortization, not merely a claim that a particular feature set predicts one
benchmark well.

## Hostile-reader checks

| Objection | Revised answer | Residual risk |
|---|---|---|
| “Predicted physical descriptors already exist.” | The Introduction cites ml-QM-GNN and the direct 2026 two-stage hydration precedent before stating the SolvAI question. | Low. |
| “Solvent-conditioned descriptors already exist.” | ML-PCM, 3D-RISM and ImPerHam are cited and distinguished by their per-query physical computation. | Low. |
| “This is feature engineering.” | Matched removal and shuffled-prior controls isolate molecule-aligned response information; block ablations show that the full response vocabulary, not any arbitrary extra columns, matters. | Moderate; the conceptual value depends on accepting these controls as evidence of an intermediate representation. |
| “PIMD was distilled.” | Introduction, figure legend, abstract and cover letter state that PIMD8 is a comparator and no PIMD feature is retained. | Low. |
| “The result is universal.” | The paper reports five-partition performance and the weaker global family/scaffold results; it restricts the domain to neutral small-molecule hydration. | Moderate but appropriately disclosed. |
| “The novelty is only these 15 features.” | The manuscript frames the contribution as a controlled response-layer experiment and a transfer criterion: relevance, structure-learnability and complementarity. | Moderate; this is the correct largest supported claim. |

## Final wording boundary

The strongest supported statement is:

> A heterogeneous layer of structure-predicted, property-proximal solvent responses
> supplies complementary information to an experimentally supervised hydration model,
> reaches the PIMD8 accuracy scale on the ARROW reference chemistry without per-query
> physical calculation, and fails when the desired response cannot be inferred
> accurately enough from structure.

The manuscript should not claim that SolvAI invented surrogate-predicted physical
descriptors, first combined physics and ML for hydration, or distilled PIMD into the
selected model.
