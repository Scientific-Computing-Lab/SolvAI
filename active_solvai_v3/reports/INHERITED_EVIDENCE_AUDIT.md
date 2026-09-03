# Active SolvAI v3 inherited-evidence audit

**Audit date:** 2026-09-03  
**Workstream:** Active SolvAI v3: Adaptive Simulation Effort  
**Branch:** `active-solvai-v3-effort-allocation`  
**Parent commit:** `770e13ce68ca80104557fd67a224abc5a2c44767`

This audit was completed before fitting or scoring any v3 difficulty model and
before launching any new molecular simulation. The earlier results are treated
as immutable evidence, including their negative conclusions.

## 1. Authoritative specification and immutable history

The operating document is
`Active_SolvAI_v3_Adaptive_Simulation_Effort_Blueprint.docx`, SHA-256
`264f7f4afad7a9fa578992b33e6d5fefa83fe15f74be7fb9d85fd7e2f751b4de`.
The complete document was extracted and inspected: 1,945 text lines, 14,040
words, all tables and equations, three embedded figures, references, Part XIII
and the Part XIV deliverables checklist. The ingestion record is
[`BLUEPRINT_INGESTION.md`](BLUEPRINT_INGESTION.md).

The inherited milestones are:

| Scientific state | Immutable commit | Status carried into v3 |
|---|---|---|
| Parent SolvAI release | `531f6cfd21e319c951b461c9ef24fa754790f91d` | Structure-only hydration model and confirmatory controls |
| Active SolvAI v1 no-go | `8fb984c2eb26d016c6b81cf488f88dc667ca9cd3` | Actual PIMD2 response did not improve the endpoint; the placement policy lost to uniform sampling |
| Active SolvAI v2 independent-noise audit | `f2cc60fb73416e0e417c97c600d5abd00868031c` | The apparent same-curve oracle headroom reversed under complementary-block evaluation |
| Active SolvAI v2 independent-replica gate | `770e13ce68ca80104557fd67a224abc5a2c44767` | Proposed eight-molecule replica experiment had effectively zero power and was not launched |

V3 does not reinterpret these findings. It changes the estimand from choosing
lambda locations to allocating additional time over a fixed conservative grid.

## 2. Canonical numerical reproduction

`scripts/reproduce_inherited.py` recomputed the values below from molecule-level
or trajectory-level machine-readable outputs; it did not copy report prose.
The complete table and hashes are in
`results/inherited/canonical_reproduction.csv` and `.json`.

| Evidence | Reproduced result |
|---|---:|
| Parent matched structure-only ARROW MAE | 0.303348 kcal mol-1 |
| Parent full SolvAI ARROW MAE | 0.202234 kcal mol-1 |
| Structure-only, five repeat mean +/- s.d. | 0.312953 +/- 0.004207 kcal mol-1 |
| Full SolvAI, five repeat mean +/- s.d. | 0.207373 +/- 0.004443 kcal mol-1 |
| Zero-ARROW-label structure-only / SolvAI MAE | 0.384822 / 0.256941 kcal mol-1 |
| Tier-A endpoint-disjoint, n=220, structure-only / SolvAI MAE | 1.531647 / 1.152554 kcal mol-1 |
| Tier-A strict response-source-disjoint, n=97, structure-only / SolvAI MAE | 2.138304 / 1.535595 kcal mol-1 |
| V1 frozen SolvAI / actual-PIMD2-residual endpoint MAE | 0.186374 / 0.189856 kcal mol-1 |
| V1 Active BQ / uniform-direct integral MAE, 5 windows | 1.701248 / 1.152880 kcal mol-1 |
| V1 same-curve development oracle integral MAE, 5 windows | 0.336580 kcal mol-1 |
| V1 Active BQ / uniform-direct integral MAE, 7 windows | 1.607515 / 1.091747 kcal mol-1 |
| V1 same-curve development oracle integral MAE, 7 windows | 0.068266 kcal mol-1 |
| V2 cross-fitted oracle / uniform-direct MAE, 5 windows | 1.732079 / 1.301076 kcal mol-1 |
| V2 cross-fitted oracle / uniform-direct MAE, 7 windows | 1.681244 / 1.272475 kcal mol-1 |
| V2 proposed-replica power, 5 / 7 windows | 0.00000 / 0.00001 |
| V2 dense-reference reliability probability | 0.44602 |

The v1 actual-response comparison worsened by 0.003482 kcal mol-1 rather than
improving the experimental endpoint. At the integral level, the v1 Active
policy was worse than uniform direct integration by 0.548369 and 0.515768 kcal
mol-1 at five and seven windows. The v2 complementary-block oracle was also
worse than uniform by 0.431004 and 0.408769 kcal mol-1. These are binding
negative results, not optimization targets for v3.

## 3. Trajectory census

`scripts/inventory_trajectories.py` searched the shared workspace, classified
every `SYSTEM` energy file, joined historical attempts to their recorded QC
status and recorded file hashes. The row-level census is in
`data/manifests/trajectory_inventory.csv` and `.parquet`.

| Inventory item | Count or status |
|---|---:|
| All located `SYSTEM` energy files | 418 |
| Qualified PIMD2 files | 363 |
| Historical probe attempts / successful files | 228 / 219 |
| Historical unique molecules / complete three-point molecules | 76 / 72 |
| V1 dense full-grid molecules | 12 |
| V1 dense manifest rows | 180 = 12 molecules x 15 lambda values |
| Reused historical windows / newly run windows | 36 / 144 |
| Production per dense molecule-window | 5 ps |
| Independent streams per dense molecule-window | **1** |
| Qualified independent replica sets | **0** |

The 15 lambda values are 0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 and 1.00. The dense library
contains Propane, Ethanol, Acetone, Pyridine, Cyclohexane, Octane,
Octanol, AceticAcid, Acetamide, DiMethylEther, EthylAcetate and
Anthracene. Each qualified 5-ps trajectory contains 50 post-initial frames at
0.1-ps spacing. The historical inventory includes one incomplete failed output;
it remains in the audit table and is not marked scientifically eligible.

Other located files do not repair the independence deficit: 21 legacy outputs
lack a qualifying compatible configuration, 22 files are repository fixtures,
10 are single-window RTX 3090 protocol benchmarks and one is ancillary. The
parent table contains 125 scalar hydration rows (85 water and 40 cyclohexane)
with classical/PIMD4/PIMD8 endpoint fields, not per-window time series.

## 4. Exposure and permitted use

The following are permanently **development-exposed** for v3:

- all 85 ARROW reference molecules and their endpoint labels;
- all Tier-A molecules and the Tier-A evaluation results;
- the 76 historical PIMD2-probe molecules, including the 72 complete
  three-lambda cases;
- the four v1 calibration and eight v1 prospective dense sentinels (the full
  twelve-molecule dense library);
- every curve, block, policy prediction and aggregate reported in v1 or v2.

These data may support provenance reproduction, diagnostic development,
power calculations and policy training under molecule-grouped validation.
They cannot constitute a new prospective sentinel or an independent external
test. Contiguous blocks or complementary halves from the single inherited
trajectory may be used only as explicitly labelled development diagnostics;
they are not independent replicas. Experimental hydration labels are barred
from v3 policy inputs, molecule selection, stopping, hyperparameter selection
and primary evaluation.

No inherited dataset can by itself validate a deployable effort-allocation
policy: no molecule-window has two genuinely independent compatible full-grid
streams. A new campaign is therefore permissible only if the frozen Gate-1
identifiability, reference-reliability and prospective-power analyses quantify
a sufficient design and its cost.

## 5. Artifact integrity

| Artifact | SHA-256 |
|---|---|
| Parent paper PDF | `7bf883c62828e08f25b65fded21a4a149dc776f198f0a26ee03394b6fd52f041` |
| Parent Supplement PDF | `58b46b69f5b6d23a21e0eece00103b9c904f4227874aad345a6bcdc87cc81d89` |
| Final SolvAI model manifest | `d64d6cbca36b9b5cc672c9b7622c5b72dc1e9409cfd0daaf6a7a6d36c5ec64c7` |
| V1 phase-1 canonical metrics | `aef996e0a1785af037221d83204d6b18644249f66d9b7647e5b294f36f78a3bb` |
| V1 dense replay canonical metrics | `e2f3558888527d51098ff89b9d3b0206da38fd7f23c183013e23e4fabc2d4717` |
| V2 independent-noise canonical metrics | `c0f01f65f3cfae4230d01fb536b7f375a1a0db29742a5023870ed53022cec2f4` |
| V2 power analysis | `70a486901e85eba0828d8a94f479c702c860acdebcf131e033b6177259de5381` |
| Tier-A predictions | `aa1caca2781137bec2b25369df93e70b56eb26cb91c10c009265e43c62139ad8` |
| V3 reproduced metrics | `c823fe726f1a1ffeb053f02201fe801bf03b91de9020b96e4baf11e8ce7f8ede` |
| V3 trajectory inventory CSV | `a052283799354b8ef4096f370944fc74d1343ae6a3c678d4030dbe8ed1ee3065` |
| V3 trajectory inventory summary | `46e4391ea12d41da90b63240faf0c84cefb5e7bf0e3eda0ae0b7fef407d44ec2` |

## 6. Environment and resource constraint

The host provides one NVIDIA RTX 3090 (24,576 MiB, driver 535.309.01), an
Intel i7-4930K with six physical/twelve logical cores, 62 GiB RAM and 116 GiB
free on a 1.8-TiB filesystem at audit time. Disk use is already 94%; any later
simulation freeze must include a storage forecast and cleanup-free retention
plan. No v3 GPU simulation has been launched.

## 7. Gate-0 conclusion

Provenance and numerical reproduction pass. The inherited evidence is adequate
for a zero-new-simulation measurement study of early diagnostics, but not for
independent validation. Gate 1 must answer whether early diagnostics and frozen
SolvAI responses predict *complementary later-window difficulty* under strict
molecule-held-out analysis. Any positive result remains development evidence;
any new physics requires a separate power-qualified campaign freeze.
