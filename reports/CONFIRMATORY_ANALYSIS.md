# SolvAI Phase 1 confirmatory analysis

## Executive result

The preregistered matched comparison supports the narrow scientific claim that
molecule-aligned, structure-predicted response priors add useful information to the
frozen endpoint pipeline. After removing standardized-equivalent benchmark records
from affected teacher sources and refitting those teachers with original split
membership preserved, the matched descriptor-only endpoint has an MAE of
**0.303 kcal mol⁻¹**, whereas the full 15-prior SolvAI stack has an
MAE of **0.202 kcal mol⁻¹**. The paired difference is
-0.101 [-0.215, -0.020] kcal mol⁻¹ (candidate minus baseline;
95% molecule bootstrap interval).

This is a stronger causal control than the prior campaign comparison. It also
changes the publication record in two important ways: the historical 0.238606 value
was not a matched descriptor-only endpoint, and the historical 0.197047 model used
teachers filtered only by exact connectivity. The standardized-exclusion
confirmatory point estimate, **0.202**, is the appropriate primary
fixed-partition value. It is on the ARROW/PIMD8 accuracy scale
(0.205) but is not evidence of superiority.

## Predeclared endpoint controls

| Fixed feature set | MAE | RMSE |
| --- | ---: | ---: |
| A_structure_only | 0.303 | 0.586 |
| B_empirical_residual | 0.243 | 0.522 |
| C_computation_core | 0.253 | 0.422 |
| D_smd_water | 0.253 | 0.395 |
| E_confsolv | 0.299 | 0.562 |
| F_full_solvai | 0.202 | 0.318 |
| G_narrow_reference | 0.212 | 0.346 |
| H_narrow_smd_reference | 0.207 | 0.325 |

The blocks are scientific ablations, not an additive attribution. In isolation,
the empirical/residual block improves the matched endpoint; SMD alone is positive
under the original exact-connectivity teachers but becomes statistically neutral
after the more conservative standardized exclusion. ConfSolv alone is neutral. The
full block is positive, consistent with complementary information across response
coordinates rather than one universally dominant teacher.

Across the five preregistered partitions, the descriptor-only endpoint is
**0.313 ± 0.004** and
full SolvAI is **0.207 ± 0.004**
(mean ± sample s.d.). The molecule-averaged paired improvement is
-0.106 [-0.216, -0.027] kcal mol⁻¹. Full SolvAI improves on the matched
endpoint in all five partitions.

## Shuffled-prior control

Meaningful molecule--response alignment is necessary. On the fixed partition,
the mean of five shuffled-prior controls is **0.306**
MAE versus **0.202** for aligned priors; the
paired difference is -0.104 [-0.214, -0.025] kcal mol⁻¹. Shuffled priors
remain near the matched structure-only endpoint rather than reproducing the aligned
gain. This rejects the explanation that arbitrary extra columns or model width are
sufficient.

## Global chemical separation

Chemical separation was applied to every endpoint-supervised molecule, including
the 1,280-label external pool.

| Regime | Structure-only MAE | Full SolvAI MAE | Paired difference [95% CI] |
| --- | ---: | ---: | ---: |
| global_butina_0_70 | 0.349 | 0.216 | -0.134 [-0.241, -0.056] |
| global_family | 1.262 | 0.468 | -0.794 [-0.963, -0.633] |
| global_nn_0.50 | 0.359 | 0.219 | -0.140 [-0.258, -0.052] |
| global_nn_0.60 | 0.344 | 0.211 | -0.133 [-0.243, -0.050] |
| global_nn_0.70 | 0.343 | 0.211 | -0.132 [-0.239, -0.052] |
| global_nn_0.80 | 0.339 | 0.210 | -0.130 [-0.238, -0.050] |
| global_scaffold | 0.728 | 0.376 | -0.352 [-0.470, -0.242] |

The prior advantage survives every predeclared separation regime, with paired
intervals below zero. Absolute error rises sharply for leave-family and scaffold
extrapolation, so these results support transfer of the response-prior advantage,
not broad high-accuracy extrapolation.

## Zero-ARROW-label transfer

With no ARROW experimental label used in endpoint fitting, the matched structure
model gives **0.385** MAE and full SolvAI gives
**0.257**. The paired difference is
-0.128 [-0.246, -0.037] kcal mol⁻¹. This supports representation value without
fold-local adaptation, but it is not an independent external validation because
the method was developed in the context of the 85-solute set.

## Identity and chemical-distance audit

The endpoint pool contains no exact, fragment-parent, uncharged-parent or canonical-
tautomer match to the 85 reference connectivities. The expanded audit found
standardized equivalents in three teacher sources: 2 CombiSolv-QM rows, 32 MolSolv
records and 22 ConfSolv rows. Those records were excluded, the affected teachers
were refitted while preserving every remaining molecule's original source split,
and all affected endpoint analyses were repeated. Morgan-similarity 1.0 collisions
between non-identical structures are reported as fingerprint collisions, not hidden
identity matches. Full details and every pair appear under `audits/confirmatory/`.

## Claims that survived

- The 15 aligned response priors materially improve an otherwise matched endpoint.
- Arbitrarily shuffled priors do not reproduce the gain.
- The advantage survives family, scaffold, cluster and nearest-neighbour separation
  applied globally to all endpoint labels.
- The advantage remains when no ARROW label is used for endpoint fitting.
- Structure-only SolvAI reaches the ARROW/PIMD8 accuracy scale on the fixed reference
  partition, with no simulation at inference.

## Claims that weakened or died

- The historical 0.238606 number cannot be described as the matched no-prior
  baseline; the correct matched value is 0.303.
- The historical 0.197047 estimate is not the conservative publication headline
  after the expanded identity audit; the corrected fixed result is 0.202.
- Robust sub-0.20 performance is not supported. Corrected repeated performance
  centres at 0.207 kcal mol⁻¹.
- Chemical extrapolation is not solved: global family and scaffold MAEs are
  0.468 and
  0.376, respectively.
- PIMD supervision was not retained; the result concerns distillation of diverse
  solvation-response calculations, with PIMD8 serving only as the accuracy
  comparator.

## Nature Communications thesis

The confirmatory evidence supports a focused thesis: external physical calculations
can define molecule-aligned response coordinates that are learned from structure and
provide a reproducible endpoint advantage beyond the same experimental labels,
representation, model class, folds and seeds. It does not establish a universal
principle across properties or chemistry. The manuscript should therefore lead with
the controlled response-prior result and use the PIMD8 comparison as a scientifically
important scale reference, not as the causal proof.

## Publication decision and figure mapping

**Go after evidence-corrected restructuring.** The confirmatory controls support
submission to *Nature Communications* as a controlled demonstration of reusable
physical-response supervision. A reasonable pre-submission estimate is 35--55% for
editorial review and 20--35% for acceptance conditional on review; the principal
risk is breadth, not internal validity.

- Main Fig. 2: matched endpoint, molecule-level paired errors, source-block
  intervals and complete split repeats.
- Main Fig. 3: global family/scaffold/cluster separation, nearest-neighbour exclusion
  and zero-ARROW-label transfer.
- Main Fig. 4: compact versus latent response representations, PIMD2 response-head
  error and downstream consequences.

The Phase 2 rebuild should use the corrected 0.202 fixed estimate and 0.207 repeat
mean, retire the historical 0.197/0.239 contrast from headline use, distinguish
PIMD8 as comparator rather than teacher, and place complete exploratory results in
Supplementary Information and Supplementary Data.

## Execution record

All definitions are frozen in `release/CONFIRMATORY_FREEZE.md`. Valid analyses used
`uv run` with Python 3.11.15, RDKit 2026.03.5 and scikit-learn 1.7.2. Primary commands:

```bash
uv run python scripts/run_confirmatory_endpoint.py --mode primary
uv run python scripts/run_confirmatory_endpoint.py --mode repeats
uv run python scripts/run_confirmatory_endpoint.py --mode shuffle
uv run python scripts/audit_confirmatory_chemistry.py
uv run python scripts/train_standardized_exclusion_teachers.py
uv run python scripts/verify_confirmatory_teacher_refits.py
uv run python scripts/run_standardized_exclusion_endpoints.py
uv run python scripts/run_confirmatory_endpoint.py --mode shuffle --standardized-exclusion
uv run python scripts/run_confirmatory_separation.py --standardized-exclusion
uv run python scripts/summarize_confirmatory.py
```

The optional experimental-label learning curve was not run because it required 750
additional endpoint fits and was not needed to resolve any primary interpretation
rule. A cross-source teacher-fidelity regression was not run because the available
teacher targets differ in units, scales, model classes and test protocols; treating
them as commensurate would be more misleading than informative.
