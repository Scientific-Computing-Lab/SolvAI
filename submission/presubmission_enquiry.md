# Presubmission enquiry

**Proposed title:** Distilling solvation physics into AI for simulation-free
hydration free energies

High-accuracy hydration calculations resolve coupled electronic, conformational
and nuclear-quantum response through expensive molecular simulation. We ask
whether that physics can instead be learned once and reused. SolvAI trains
structure-to-response surrogates on strictly benchmark-disjoint physical
calculations, then predicts hydration free energy from SMILES alone. On the
85-solute reference set introduced with ARROW, it improves a direct
structure-only MAE from 0.239 to 0.197 kcal/mol and reaches the reconstructed
ARROW/PIMD8 accuracy of 0.205 kcal/mol. Across five independent partitions the
mean is 0.204 ± 0.005 kcal/mol. Mechanistic ablations show that aligned
water-response supervision transfers, whereas several broader physical
representations and sparse PIMD2 response curves do not. We believe the broad
interest lies in the change of computational paradigm: simulation becomes a
reusable source of physical supervision rather than a calculation required for
every prediction.
