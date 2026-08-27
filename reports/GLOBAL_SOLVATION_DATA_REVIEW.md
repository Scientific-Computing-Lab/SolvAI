# Global solvation data and model review

## Outcome of the reconnaissance

The search found **41 distinct source families** spanning experimental
hydration, million-scale solvent-conditioned QM targets, explicit alchemical
simulation results, non-covalent quantum interaction data, and pretrained molecular
encoders. The most immediately useful additions are:

1. **ConfSolv water response hierarchy:**
   5,392,567 H2O conformer records were mapped and
   aggregated into 39,878 neutral,
   benchmark-disjoint connectivity-level examples. The usable targets are the
   gas/solution conformer corrections and response-distribution summaries; raw
   COSMO-RS hydration values are not treated as calibrated experimental truth.
2. **MolSolv SMD(water):** 1,729,545 source
   conformers yielded 350,391 sampled,
   deduplicated structures for water-specific M06-2X/6-31G* SMD pretraining after
   removing 82 exact benchmark structures.
3. **Explicit-solvent phase-space averages:** 559
   strict nonbenchmark FreeSolv connectivities provide 114
   trajectory-averaged geometry/response targets for structure-only distillation.
4. **Paired gas/SMD QM–NNP hierarchy:** 330
   complete strict molecules provide solvent-induced geometry response plus QM and
   equivariant-NNP energy decompositions.
5. **CombiSolv-QM:** 1,000,000 solute-solvent pairs,
   including 3,988 water rows. The strict copy has
   992,329 pairs after globally excluding every
   ARROW-benchmark solute.
6. **G4MP2/Foundry-ML:** 130,258 molecules with five continuum
   solvation energies plus charges, dipole, polarizability, orbitals, and
   thermochemistry; 130,202 remain after benchmark exclusion.
7. **OpenFE/OpenFF ASFE:** 3,618 successful
   protocol-unit results retain 14-state MBAR overlap, forward/reverse convergence,
   and replica-exchange diagnostics. After global identity exclusion,
   523 solutes remain.
8. **DES370K water interactions:** 170,486
   neutral water–solute geometries provide coupled-cluster interaction-energy
   distributions and SAPT components for 306
   strict nonbenchmark solutes.
9. **Aligned quantum/alchemical additions:** 71
   benchmark-disjoint NQELiq structures carry paired classical/PIMD property shifts;
   32 strict Bannan solutes retain
   19 water MBAR increments each;
   and 46 strict SAMPL4 solutes carry multiple AMOEBA/GAFF protocols.
10. **OpenADMET foundation weights:** downloadable Chemprop message-passing
   checkpoints pretrained on up to 10M structures and on MiniMol/Jazzy/QM-derived
   descriptors. The MiniMol and generic MolPILE variants are clean candidates. The
   Jazzy variant is diagnostic-only until its target ancestry is isolated.
11. **Freecurve privileged labels:** classical ARROW, PIMD8, NQE corrections, 40
   cyclohexane PIMD4 values, and 72 complete short-PIMD2 three-point response
   curves. These are used only as outer-training auxiliary targets.

No public checkpoint was found that both (a) claims approximately 0.20 kcal/mol on
the same neutral-water benchmark and (b) provides enough identity/split provenance
to support a leakage-safe direct comparison. Published FreeSolv random-split scores
are not treated as equivalent to this ARROW benchmark.

## Source map

| source | molecules_or_samples | targets | physics_content | water_specific | data_class | downloaded | usable_for_training |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Freecurve ARROW hydration benchmark | 85 molecules / 125 molecule-solvent rows | experiment; classical ARROW; PIMD4; PIMD8; NQE and experimental residuals | high-fidelity alchemical free energies and bead convergence | Primary 85 water; 40 cyclohexane auxiliary | experimental + classical MD + PIMD | True | Yes; outer-fold masked labels only |
| Freecurve short PIMD2 response teachers | 72/76 at lambda=0.1; 74/76 at 0.5; 73/76 at 0.9; 72 complete curves | dH/dlambda and component/time-series/shell/bead statistics | short explicit-water PIMD2 response fingerprints and sparse alchemical curves | Yes | PIMD2 training-only teacher | True | Yes; never an inference input |
| FreeSolv v0.52 | 642 molecules | experimental and GAFF hydration free energy with uncertainty | experimental hydration plus classical alchemical estimates | Yes | experimental + MD | True | Yes after benchmark exclusion |
| CombiSolv-Exp-8780 | 8,780 solute-solvent records; 1,305 solutes | experimental solvation free energy and replicate spread | solvent-conditioned experimental response | No; 1,153 water rows | experimental | True | Yes after solute-level benchmark exclusion |
| CombiSolv-QM | 1,000,000 pairs; 11,021 solutes; 284 solvents | COSMOtherm solvation free energy | quantum-continuum solvent response across 284 solvents | No; 3,988 water rows | QM/continuum | True | Yes after strict solute exclusion |
| G4MP2 solvation/QM9 (Foundry-ML) | 130,258 molecules | five solvent energies; charges; dipole; polarizability; orbital and thermochemical properties | DFT continuum solvation plus quantum molecular observables | Includes water and four organic solvents | QM/DFT | True | Yes after strict benchmark exclusion |
| GuthrieSolv | 53,895 raw literature records; 3,316 Jazzy-curated rows | heterogeneous hydration-related measurements | experimental source records and uncertainty/context | Mostly hydration; raw table mixes properties/units | experimental | True | Curated subsets only; not used blindly |
| SAMPL1 hydration challenge | 63 molecules | experimental hydration free energy and uncertainty | blind-challenge experimental reference values and 3D structures | Yes | experimental benchmark | True | Yes after benchmark exclusion; small and out-of-domain |
| Solv@TUM | 658 solutes / 5,952 logK pairs | gas-solvent partition coefficient converted to solvation free energy | experimental multi-solvent response; dipole/polarizability | No; non-aqueous | experimental | True | Yes as solvent-conditioned auxiliary task |
| Thermodynamic-length learning dataset | 642 FreeSolv molecules / 111 relative transformations | relative hydration free energies and variances | alchemical network variance and trajectory-derived features | Yes | classical MD | True | Potentially, after pair/trajectory parsing |
| Jazzy | 292 Gerber fit molecules; 3,316 Guthrie validation rows | hydration free energy and H-bond donor/acceptor strengths | EEQ charges, atomic polarizabilities, polar/apolar decomposition | Yes | deterministic QSPR | True | Only label-free atomic strengths are clean headline features |
| SolvBERT | 1,000,000 QM + 8,780 experimental pairs | solvation free energy; solubility | solvent-solute language-model pretraining | No | pretrained architecture/data | True | Data yes; no downloadable checkpoint found |
| 3DMRL | Uses CombiSolv pretraining; checkpoint absent | molecular interaction properties | virtual 3D solute-solvent interaction environment | No | 3D pretrained architecture | True | Architecture feasible; full retraining not first-line |
| GeoMAW-Solv / Solvaformer | Code; external BigSolDBv2/CombiSolv | solubility and solvation tasks | 3D Equiformer/MPNN; optional learned charges | No | 3D/graph architecture | True | Code usable; advertised checkpoints require W&B access |
| DES370K / DES5M | 170,486 neutral water-dimer geometries / 306 strict solutes; 370,000 total DES370K geometries | CCSD(T), SNS-MP2, MP2, HF, SAPT0 interaction energies | water-interaction distributions plus SAPT electrostatic, exchange, induction, and dispersion components | Yes for the extracted teacher subset | QM dimers | True | Yes; aggregated structure-only teacher built |
| Uni-Mol molecular encoder | 209M conformers in published pretraining; 1,232 local embeddings | 3D denoising and molecular representation pretraining | geometry-aware atom and molecular representations | No | pretrained 3D foundation encoder | True | Yes: frozen embeddings from deterministic RDKit conformers |
| SPICE | Large conformer/interaction collection (version-dependent) | DFT energies, forces, charges, dipoles, bond orders | intramolecular and noncovalent quantum response | No; includes water-containing subsets | QM conformers/interactions | False | Potential foundation pretraining; indirect for hydration |
| AqSolDB / BigSolDBv2 | 9,982 compounds / large multi-solvent solubility corpus | aqueous or solvent-conditioned logS | solubility (conflates hydration and condensed-phase effects) | AqSolDB yes; BigSolDB no | experimental solubility | False | Only weak auxiliary representation supervision |
| ChemBERTa-77M-MTR / MoLFormer-XL | 77M / approximately 1.1B pretraining molecules | masked-token molecular representation | broad chemical structure only | No | pretrained foundation encoders | True | Yes: frozen embeddings then bounded fine-tuning |
| OpenADMET Chemprop foundation weights | 1M PubChem pretraining molecules; multiple 34 MB checkpoints | Jazzy, MiniMol, Mordred, ECFP and other pseudo-descriptors | graph encoders distilled from quantum/property representations | Jazzy checkpoint is hydration-oriented | pretrained graph checkpoints | True | Yes; use as initialization with a leakage caveat for Jazzy |
| SoluteML dGsolvDB3 / SoluteDB | 5,075 strict water hydration values; 8,098 Abraham-response molecules | experimental hydration dG; Abraham E/S/A/B/L | water partition response and hydrogen-bond/cavity/polarity axes | Hydration values yes; Abraham axes transferable | experimental + empirical physical descriptors | True | Yes after strict benchmark exclusion |
| CombiSolvH-QM / SolProp v1.2 | 2,816 strict water solutes | COSMO-RS dH; matched dG; T*dS=dH-dG | water-specific thermodynamic response hierarchy | Yes for processed teacher subset | QM/continuum thermodynamics | True | Yes after strict benchmark exclusion |
| DISSOLVE2-ANIONS | 3,662 strict neutral acid structures | aqueous pKa; gas-phase acidity; COSMO-RS anion solvation | deprotonation and charged water-response hierarchy | Yes | experimental + QM/continuum | True | Safe targets only; upstream neutral dG quarantined |
| Relative-solvation multi-code archive | 9 absolute solutes plus relative transformations | lambda-resolved gradients for archived runs; multi-code dG | explicit alchemical response under AMBER/CHARMM/GROMACS/SOMD | Yes | classical explicit-solvent MD | True | Too few distinct solutes for a general lambda model |
| Thermodynamic-length relative hydration pairs | 79 strict pairs / 94 molecules | vacuum/solvent relative dG and uncertainty | phase-resolved relative alchemical response | Yes | classical explicit-solvent MD | True | Yes as a small pairwise/contrastive auxiliary task |
| SAMPL4 AMOEBA polarizable hydration protocols | 51 strict SAMPL4 structures; 46 with seven triplicate AMOEBA/GAFF protocols | AMOEBA Poltype/polarization-group/multiconformer/basis/water/OH-scale dG; GAFF dG; protocol spread and deltas | polarizable versus fixed-charge explicit-water alchemical response | Yes | classical polarizable/fixed-charge MD | True | Yes as benchmark-disjoint privileged targets |
| Bannan et al. partition-coefficient alchemical archive | 32 strict of 41 solutes; 19 water alchemical edges per solute | 19 adjacent-state MBAR increments and uncertainties; GAFF/DC totals; curve PCs; water/cyclohexane/octanol archive | explicit-solvent alchemical response curves | Processed teacher is water; archive includes two organic solvents | classical explicit-solvent MD | True | Yes as a small lambda-response teacher |
| AMOEBA multi-solvent small-molecule solvation | 21 solutes in toluene/chloroform; 6 in acetonitrile/DMSO | experimental, triplicate AMOEBA, and triplicate GAFF solvation dG | polarizable/fixed-charge response across four nonaqueous solvents | No | experimental + classical polarizable/fixed-charge MD | True | Potential solvent-conditioned auxiliary task; too small alone |
| Replica Exchange with Flexible Timing hydration benchmark | 7 molecules / 112 convergence rows | hydration dG and uncertainty versus method and ns/window | sampling/convergence response across REFT, FEP, and Transformato | Yes | classical explicit-water alchemical MD | True | Diagnostic convergence supervision; only seven molecules |
| OpenFE/OpenFF 2.3.0 FreeSolv ASFE archive | 603 solutes x 3 repeats x solvent/vacuum legs; 523 strict nonbenchmark solutes | hydration and leg dG; uncertainties; 14-state MBAR overlap matrices; forward/reverse convergence; replica-exchange spectra; NAGL charges | 10 ns explicit-water alchemical response, convergence, and mixing diagnostics | Yes | classical explicit-water alchemical MD | True | Yes as benchmark-disjoint privileged supervision |
| NQELiq-298 | 92 molecular liquids; 71 strict nonbenchmark structures | classical/PIMD-H/PIMD-D density, volume, expansivity, compressibility, dielectric constant, dHvap; NQE and isotope shifts | paired classical and path-integral bulk-liquid response | No; chemically diverse neat molecular liquids | classical MD + PIMD | True | Yes as benchmark-disjoint NQE auxiliary supervision |
| MolSolv SMD(water) | 1,729,545 source conformers; 350,391 strict sampled structures | M06-2X/6-31G* SMD(water) solvation free energy | large water-specific quantum-continuum solvent-response teacher | Yes | DFT/SMD continuum water | True | Yes after global benchmark exclusion |
| FreeSolv explicit-solvent phase-space averages | 642 source molecules; 559 strict nonbenchmark connectivities | Boltzmann-averaged intramolecular distance matrices | explicit-water trajectory-averaged geometry, conformational response, distance/Coulomb spectra | Yes | explicit-water MD phase-space averages | True | Yes as benchmark-disjoint privileged supervision |
| G-NequIP paired gas/SMD solvation benchmark | 428 energy records / 424 paired geometries; 355 strict nonbenchmark connectivities | QM SMD dG; NNP solvent/geometry/total response; solvent-induced geometry changes | paired gas-to-water implicit-solvent energy hierarchy | Yes | DFT/SMD + equivariant NNP | True | Yes as benchmark-disjoint privileged supervision |
| MLFF hydration free-energy comparison | 59 source solutes; 45 strict rows | Organic-MPNICE, GAFF, OPLS4, E-sol, and DFT/PBF hydration dG | matched multi-fidelity force-field hydration hierarchy | Yes | MLFF + classical force fields + continuum DFT | True | Yes, but only 45 strict structures |
| FreeSolv archived per-window GROMACS energies | Documented for the 2017 642-molecule campaign | 20-state fep/vdW dH/dlambda and XVG energy time series | full classical alchemical response curves | Yes | classical explicit-water alchemical MD | False | Not currently: documented directory absent from archives |
| GNNImplicitSolvent explicit-water mean-force model | 369,486 molecules / approximately 3.2M conformers; 5,156 local static evaluations | explicit-water mean solvation forces and learned solvent energy | force-matched solvent response from explicit-water trajectories | Yes | explicit-water MD force teacher + public GNN checkpoint | True | Yes; checkpoint/static response is benchmark-label independent |
| Lambda-aware implicit-solvent force checkpoint | 280K training configurations reported; 1,363 local static evaluations | solvent forces plus steric/electrostatic lambda derivatives | learned lambda-resolved solvent response and endpoint energy | Yes | explicit-solvent force/derivative distillation checkpoint | True | Yes; deterministic structure/conformer evaluation |
| ConfSolv COSMO-RS conformer solvation | 5,392,567 H2O conformer rows; 39,878 strict neutral connectivities | Boltzmann-ensemble hydration plus conformer-resolved COSMO-RS response | gas/solution conformer corrections and water response distributions | Yes (H2O slice extracted from 41-solvent archive) | DFT/COSMO-RS conformer ensemble | True | Yes; structure mapping and benchmark exclusion complete |
| ReSolv implicit-solvent free-energy potential | 583 paired vacuum/water trajectories; 389 FreeSolv training molecules | experimental hydration free energy through differentiable BAR training | vacuum/water conformer distributions and learned solvent potential | Yes | ML potential + trajectory reweighting | True | No for headline benchmark; checkpoint is target-contaminated |
| Freecurve all-branch Git physics audit | 89 branches / 7,513 physics-related path matches | TI/BAR/dHdl/PIMD/energy output discovery | historical Git-tree audit across priority repositories | Mixed | repository provenance audit | True | Only identified molecule-resolved outputs |

## Search coverage and retrieval

The scan covered Crossref/OpenAlex-indexed literature, GitHub code/repository search,
Zenodo/DataCite metadata, Hugging Face model and dataset indices, publisher
supplements, FreeSolv/GuthrieSolv/SAMPL archives, and the available Freecurve Git
history. Query responses are preserved in `data/catalog/raw/`; downloaded repositories
and model cards are under `data/external/`.

High-value artifacts inspected recursively include SolvBERT/CombiSolv, 3DMRL,
GeoMAW-Solv/Solvaformer, Hydra, thermodynamic-length learning, qmmm-hydration,
Jazzy, OpenADMET foundation-model code and weights, DES370K/DES5M, SPICE,
Solv@TUM, SAMPL1, AqSolDB/BigSolDB, ChemBERTa, MoLFormer, Uni-Mol, the
Weinreich phase-space archive, and the G-NequIP paired gas/SMD benchmark.

The FreeSolv 2017 README and its 20-state GROMACS protocol were also traced to
GitHub, Zenodo, eScholarship, and the maintainer's issue response. The documented
`gromacs_energies` XVG directory is absent from the downloadable GitHub/Zenodo ZIP,
so the full 642-molecule lambda curves remain a specifically identified data gap.

## Scientific triage

- **Use now:** strict CombiSolv-QM water and solvent-conditioned pairs; strict
  G4MP2 multi-property targets; clean public hydration labels; MiniMol/MolPILE
  pretrained graph encoders; label-free Jazzy atomic strengths; ARROW/PIMD targets.
- **Use now:** the benchmark-excluded MolSolv SMD(water) corpus is the largest
  water-specific quantum-continuum teacher found and is used for structure-encoder
  pretraining, not as an inference-time calculation.
- **Use now:** ConfSolv contributes a much larger water-specific conformer-response
  hierarchy. Its conformer corrections and distribution summaries are distilled
  through structure-only surrogates; its poorly calibrated absolute hydration
  target is not mixed indiscriminately with experiment.
- **Use now:** explicit-water phase-space averages and paired gas/SMD QM–NNP
  differences provide narrow dynamic-geometry and solvent-response targets for
  structure-only students.
- **Use as privileged supervision only:** PIMD2/PIMD4/PIMD8, NQE residuals,
  component energies, alchemical statistics, and classical ARROW labels.
- **Diagnostic only:** Jazzy's fitted hydration-energy outputs and its public
  foundation checkpoint, because 76/85 benchmark connectivities occur in the
  292-molecule Gerber coefficient-fit table.
- **Use now:** the extracted DES370K water subset provides compact per-solute
  coupled-cluster/SAPT response targets without loading its geometry archive during
  student training. DES5M and full SPICE retraining remain higher-cost follow-ups.
- **Use now:** OpenFE ASFE diagnostics, Bannan response curves, SAMPL4 AMOEBA
  protocols, and NQELiq classical/PIMD shifts are benchmark-disjoint, aligned
  physical teachers. Their small or domain-shifted source sizes are explicitly
  tested rather than assumed to transfer.
- **Not blindly combined:** MNSol/Guthrie and multi-solvent experimental tables mix
  methods, standard states, temperatures, and measurement types. They require
  solvent/method conditioning and identity exclusion.

## Key references

- FreeSolv: DOI `10.1021/ci3001277`, update `10.1021/acs.jced.7b00104`.
- ARROW/PIMD benchmark: DOI `10.1038/s41467-022-28041-0`.
- CombiSolv transfer learning: DOI `10.1016/j.cej.2021.129307`.
- ConfSolv conformer ensembles: DOI `10.5281/zenodo.8292520`.
- MolSolv SMD(water): DOI `10.5281/zenodo.7262826`.
- G4MP2 solvation collection: DOI `10.18126/jos5-wj65`.
- Jazzy: DOI `10.1038/s41598-023-30089-x`.
- DES370K/DES5M: DOI `10.1038/s41597-021-00833-x`, Zenodo 5676266/5706002.
- SPICE: DOI `10.1038/s41597-022-01882-6`.
- SAMPL1 hydration challenge: DOI `10.1021/jp806724u`.
- Hybrid alchemical ML: DOI `10.1021/acs.jcim.0c00600`.
- Learning remaining force-field error: DOI `10.1021/acs.jctc.3c00981`.
- Kernel hydration prediction and database bias: DOI `10.1063/5.0012230`.
- Physically inspired hydration descriptors: DOI `10.1021/acs.jpclett.2c03858`.

## Important limitations

“Downloaded” does not imply unrestricted redistribution: each row records the
available license or publisher terms. Data with unclear licensing remain local and
are not bundled into a deployable artifact. The public scan is broad but cannot
prove nonexistence of unindexed/private data. Dataset sizes and methods are not
accuracy-equivalent; only predictions produced on the fixed 85-molecule benchmark
under its stored folds enter the headline comparison.
