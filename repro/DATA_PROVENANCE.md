# Data provenance

All external supervised records were canonicalized before fitting. The confirmatory
release excludes exact connectivity matches and standardized equivalents obtained by
cleanup, largest-fragment selection, uncharging and canonical tautomerization. Exact
audit records are under `audits/confirmatory/`; historical source hashes remain in
`data/manifests/training_source_manifest.json`. Counts below distinguish original
source units from the structures used by the final response teachers.

| Source | Scientific role | Original size | Filtered size | Benchmark overlaps removed | Version | URL / DOI | Processed SHA-256 | Licence | Redistributed? | Reproduction command |
|---|---|---:|---:|---:|---|---|---|---|---:|---|
| MolSolv SMD(water) | water-specific quantum-continuum response | 1,729,545 calculations; 710 MB archive | 350,359 structures | 32 standardized equivalents beyond the original exact exclusion | Zenodo 7262826 | [10.5281/zenodo.7262826](https://doi.org/10.5281/zenodo.7262826) | `758b9b56a991d874fb1067eca024dbc2244e67f8c8eb319db08ec4ab33e07460` (historical processed table; confirmatory exclusions separately hashed) | CC BY 4.0 data | yes | `python scripts/download_data.py molsolv` |
| ConfSolv H2O | conformer and solvent-response hierarchy | 5,392,567 H2O conformer records; 11.5 GB archive | 17,829 complete teacher structures (39,878 broader filtered connectivities) | 22 standardized equivalents beyond the original exact exclusion | Zenodo 8292520 | [10.5281/zenodo.8292520](https://doi.org/10.5281/zenodo.8292520) | `d283d2da04987a1802dbf60c5606394ac745c26609e2628144874b27bd037e27` (historical processed table; confirmatory exclusions separately hashed) | CC BY 4.0 | yes | `python scripts/download_data.py confsolv` |
| SoluteML | experimental hydration and Abraham response axes | 9.1 MB data archive | 3,928 added hydration structures and 8,098 Abraham structures | 84 connectivities | Zenodo 5792296 | [10.1021/acs.jcim.1c01103](https://doi.org/10.1021/acs.jcim.1c01103) | `4b90805c619a53b3b25a479c5c6b393333b575473738d5b290d899fb4a3e5fa2` (Abraham table) | CC BY 4.0 | yes | `python scripts/download_data.py soluteml` |
| OpenFE/OpenFF 2.3.0 FreeSolv ASFE | explicit-water alchemical response | 603 solutes; 7.9 MB archive | 520 structures | 83 | Zenodo 21810272 | [10.5281/zenodo.21810272](https://doi.org/10.5281/zenodo.21810272) | `3d34e27040d69fe01c45b0b9e6424baee2ba9e3c220bc6659a89d30232f9b087` | CC BY 4.0 | yes | `python scripts/download_data.py openff` |
| GBn2 / GNNImplicitSolvent | implicit and learned explicit-water response | 550 source structures used | 550 | 0 | catalogued repository version | [D4SC02432J](https://doi.org/10.1039/D4SC02432J) | `1a50219f3a0c0c60a754b430b4bc9c6c26beea1ed198e76e8aef7b14a8fa13a4` | CC BY-SA 4.0 data | yes | `python scripts/download_data.py implicit` |
| CombiSolv-QM water | COSMOtherm water response | 3,988 water records | 3,959 structures | 2 standardized equivalents beyond the original exact exclusion | publisher supplementary file | [10.1016/j.cej.2021.129307](https://doi.org/10.1016/j.cej.2021.129307) | `6fbaceff441b62c8399aae2dbd3b6cfbdfa10b86fc47cee8d60001c73da24832` (historical processed table; confirmatory exclusions separately hashed) | publisher supplementary terms; no standalone data licence | no | `python scripts/download_data.py combisolv-qm --accept-source-terms` |
| Public hydration endpoint labels | experimental hydration-response teacher | 5,075 merged structures | 1,280 selected labels | benchmark connectivities removed before merge | frozen source versions in endpoint manifest | FreeSolv, CombiSolv-Exp and SoluteML (see manifest) | `603ed02b6be25d9a3057e321f2c6ea135b012666cfdb8a1b160e37f347951ec4` (expanded table) | mixed; includes publisher-supplement records without standalone licence | no | `python scripts/download_data.py soluteml` then follow `FULL_REPRODUCTION.md` |

## Benchmark and predictions

`data/benchmark/arrow_solvation_master.parquet` is the canonical 85-solute
analysis table. Its source columns trace values to recovered Freecurve files and
the ARROW publication. The table is included because every manuscript result
depends on it and the molecular values already appear in the associated
publication materials. `results/predictions` contains held-out model outputs,
not deployment-refit predictions.

## Redistribution policy

The two large CC BY 4.0 response tables included here are processed,
benchmark-disjoint derivatives whose size is modest enough for Git LFS. Original
multi-gigabyte archives remain at Zenodo. CombiSolv-QM and the merged public
hydration tables are not copied because publisher-supplement components expose
no standalone data licence. Their trained response checkpoints, exact processed
hashes, counts, source links and reconstruction routes are preserved.

No private Freecurve repository, credential, full trajectory, raw PIMD output or
third-party archive is included in the public release.
