# Novelty positioning audit

## Question

What is the largest defensible SolvAI contribution after giving the closest prior
work full credit?

## Closest conceptual precedents

| Prior line | Representative primary work | What was already established | Relation to SolvAI |
|---|---|---|---|
| Direct hydration learning | Rauer and Bereau, *J. Chem. Phys.* (2020), [doi:10.1063/5.0012230](https://doi.org/10.1063/5.0012230) | Molecular representations can learn hydration endpoints; database composition and chemical coverage matter. | SolvAI is not novel because it predicts hydration from structure. |
| Computed-to-experimental transfer | Vermeire and Green, *Chem. Eng. J.* (2021), [doi:10.1016/j.cej.2021.129307](https://doi.org/10.1016/j.cej.2021.129307) | Large calculated solvation data can improve an experimental free-energy model. | SolvAI is not the first use of computed supervision for experimental solvation. |
| Single-calculator solvation surrogate | Ward *et al.*, *J. Phys. Chem. A* (2021), [doi:10.1021/acs.jpca.1c01960](https://doi.org/10.1021/acs.jpca.1c01960) | A graph model can replace an implicit-solvent calculation and predict its outputs from SMILES. | SolvAI does not merely emulate one calculated solvation endpoint; it uses several predicted response coordinates in a separately experimentally supervised model. |
| Physical descriptors for hydration | Zhang *et al.*, *J. Phys. Chem. Lett.* (2023), [doi:10.1021/acs.jpclett.2c03858](https://doi.org/10.1021/acs.jpclett.2c03858) | Electrostatic, polarizability, surface and related descriptors can improve hydration prediction. | SolvAI is not the first physics-informed hydration model. |
| Physical representation computed per query | Alibakhshi and Hartke, *Nature Communications* (2021), [doi:10.1038/s41467-021-23724-6](https://doi.org/10.1038/s41467-021-23724-6); Subramanian *et al.*, *JCIM* (2020), [doi:10.1021/acs.jcim.0c00065](https://doi.org/10.1021/acs.jcim.0c00065); Alibakhshi and Hartke, *Nature Communications* (2022), [doi:10.1038/s41467-022-28912-6](https://doi.org/10.1038/s41467-022-28912-6) | Continuum components, 3D-RISM hydration thermodynamics and solvent-perturbed Hamiltonian attributes are informative ML representations. | These methods establish the value of solvent-conditioned response, but retain a physical calculation for the query molecule. SolvAI predicts its response layer from structure. |
| Simulation-assisted correction | Scheen *et al.*, *JCIM* (2020), [doi:10.1021/acs.jcim.0c00600](https://doi.org/10.1021/acs.jcim.0c00600) | ML can correct alchemical free-energy calculations. | SolvAI does not require the alchemical result at inference. |
| Surrogate-predicted physical descriptors | Stuyver and Coley, *J. Chem. Phys.* (2022), [doi:10.1063/5.0079574](https://doi.org/10.1063/5.0079574) | A first model can predict QM descriptors that are passed with structure to a downstream predictor. | The two-stage architecture is established and must not be claimed as novel. |
| When predicted physics helps | Li *et al.*, *JACS* (2024), [doi:10.1021/jacs.4c04670](https://doi.org/10.1021/jacs.4c04670) | QM descriptors help chiefly when they correlate with the endpoint and can be computed or inferred accurately. | SolvAI confirms these conditions for solvent response and adds a directly tested complementarity requirement. |
| Descriptor versus surrogate latent | Chen and Stuyver, *Digital Discovery* (2025), [doi:10.1039/D5DD00256G](https://doi.org/10.1039/D5DD00256G) | Hidden surrogate representations can outperform explicit descriptor predictions in small-data chemistry tasks. | SolvAI finds the opposite ordering for ConfSolv: compact response summaries beat graph and feed-forward latent coordinates. The combined literature supports a task-dependent, not universal, choice. |
| Two-stage hydration prediction | Jia *et al.*, *Chemistry--Methods* (2026), [doi:10.1002/cmtd.202500150](https://doi.org/10.1002/cmtd.202500150) | Thirteen physically inspired solute descriptors can be predicted from a graph and supplied to a hydration model, removing their electronic-structure calculation. | This is the closest architectural precedent. SolvAI must not claim priority for predicting physical descriptors before a hydration endpoint. |

## Defensible novelty

The novelty does **not** reside in direct hydration regression, using physical
descriptors, computed-to-experimental transfer, or a generic two-stage
descriptor-surrogate architecture. The supported contribution is:

> SolvAI constructs a heterogeneous, property-proximal solvent-response layer from
> distinct solvation formalisms and demonstrates, with matched and destructive
> controls, that this structure-predicted layer carries complementary information into
> an experimentally supervised hydration model after all source calculations have
> been removed from inference.

Three parts of that statement are experimentally isolated:

1. **Response, not capacity:** the matched no-prior endpoint changes only the 15
   response columns; shuffled priors abolish the gain.
2. **Transfer, not local endpoint adaptation:** the advantage persists across five
   partitions, global family/scaffold/cluster/similarity separation and zero-ARROW-label
   transfer.
3. **A bounded mechanism:** high physical fidelity is insufficient. A response must
   be endpoint-relevant, inferable from structure and complementary to the deployed
   representation. Non-improving PIMD2 and latent-response experiments define this
   boundary.

## Wording decisions

- Retain: “heterogeneous solvent-response layer”, “structure-predicted response
  priors”, “simulation-free inference”, “PIMD8-level accuracy”.
- Avoid: “third arrangement”, “first two-stage hydration model”, “distils PIMD”,
  “new paradigm”, “physics descriptors are novel”, or “physical computation is absent
  from training”.
- Treat PIMD8 only as the high-fidelity comparator. No PIMD-trained feature is present
  in the selected model.
- Present “simulation as reusable supervision” as the broader implication supported by
  this controlled solvation demonstration, not as a universal priority claim.

## Editorial conclusion

The paper needed to exist because the prior literature established the ingredients
but did not answer the controlled solvent-response question. SolvAI joins several
distinct response formalisms in one explicit intermediate layer, tests whether that
layer adds molecule-aligned information rather than feature capacity, and shows where
amortization fails when response inference is too noisy. This is a scientific result
about the transfer of solvation physics, not merely a new descriptor stack.
