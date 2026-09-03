# Data provenance

This living file distinguishes private raw data, redistributable derived data and public source data. Exact hashes and row-level exposure are written under `data/manifests/` and `data/identity/`.

| Source | Scientific role | Local location | Redistribution status | Current use |
|---|---|---|---|---|
| Frozen SolvAI release | Prior and endpoint reference | parent repository at commit `531f6cfd21e319c951b461c9ef24fa754790f91d` | Public release | Phase 0 reproduction and all matched comparisons |
| ARROW-85 benchmark | Experimental endpoint, classical/PIMD4/PIMD8 references | `data/benchmark/arrow_solvation_master.parquet` | Released processed benchmark | Phase 1 endpoint scoring and identity mapping |
| Three PIMD2 5 ps campaigns | Actual target-molecule response at nominal lambda 0.1/0.5/0.9 | workspace `results/physics_probes/pimd2_lambda*_5ps/` | Private raw simulation output; do not redistribute | Phase 1 actual-observation gate |
| SolvAI multilambda tables | Historical predicted-PIMD2 analysis | parent `results/ablations/` and workspace `results/` | Derived tables subject to parent policy | Reproduce historical negative and map actual targets |
| Tier-A 220/97 | Development evidence only | parent `results/tier_a_external/` | Released processed tables | Diagnosis/power only; never a blind Active SolvAI test |
| Prospective dense PIMD2 sentinel | Same-Hamiltonian reconstruction and acquisition test | `simulations/dense_pimd2/` and `results/phase2/` | Raw logs/energies excluded; hashed derived response and prediction tables retained | Four calibration and eight prospective molecules on a fixed 15-point λ grid |

The original raw physical-response data remain in authorized local project storage. The
prospective release tracks configuration and energy-file hashes but does not commit raw
Arbalest logs or trajectories. Public release decisions follow the parent provenance policy
and source licences.
