# Active SolvAI final report

> **Subsequent diagnostic qualification (2026-09-03):** this report preserves
> the original prospectively registered no-go result. Its interpretation of the
> same-curve oracle as evidence of stable placement headroom was superseded by
> the independent-noise audit: held-out-block headroom reversed sign and the
> selected lambda sets were unstable. The proposed independent-replica
> resolution study was then stopped by a prospectively frozen power gate before
> simulation. See `ORACLE_INDEPENDENT_NOISE_DIAGNOSTIC.md` and
> `INDEPENDENT_REPLICA_RESOLUTION_REPORT.md`.

## Executive outcome

**No-go under the prospectively frozen program.** Actual short PIMD2 observations did not improve the experimental hydration endpoint, and the molecule-conditioned Bayesian-quadrature policy did not beat the strongest simple schedule in a prospective dense same-Hamiltonian test. The conditional multi-fidelity direction was therefore not launched.

This is a rigorous negative result for the tested 5-ps PIMD2 protocol and acquisition model. It is not evidence that all forms of active free-energy simulation are impossible.

## Direction A — experimental endpoint

On 72 molecules, the five-partition frozen SolvAI MAE was 0.186374 kcal mol⁻¹ and the actual-minus-predicted three-point residual model was 0.189856. The paired change was +0.003482 (95% CI +0.000960, +0.006093); only 31.9% of molecules improved. The aligned residual was indistinguishable from the shuffled control.

**Decision:** Direction A is killed for these observations. No λ subset, family, component or post-hoc endpoint correction was used to rescue it.

## Direction B — same-Hamiltonian reconstruction

Four calibration and eight prospective molecules were simulated on a fixed 15-point λ grid. Every molecule reused the three inherited λ=0.1, 0.5 and 0.9 windows and added 12 new 5-ps PIMD2 windows. Calibration choices were locked before any prospective response was generated.

| Windows | Active BQ MAE | 90% coverage | Strongest comparator | Comparator MAE | Paired difference | 90% CI | Molecules improved |
|---:|---:|---:|---|---:|---:|---:|---:|
| 5 | 1.701 | 0.750 | uniform_direct | 1.153 | +0.548 | [-0.101, +1.159] | 0.250 |
| 7 | 1.608 | 0.875 | uniform_direct | 1.092 | +0.516 | [-0.392, +1.491] | 0.375 |

The non-deployable same-curve oracle reached 0.337 and 0.068 kcal mol⁻¹ at five and seven windows. Subsequent independent-block scoring showed that these optimistically selected values do not establish stable molecule-specific placement; the frozen acquisition rule also did not find useful points prospectively.

**Decision:** Direction B fails the prospectively frozen criterion. Direction C was prospectively contingent on this result and was not launched.

## Failure attribution

| Possible cause | Evidence and conclusion |
|---|---|
| Lack of endpoint residual signal | Supported for this protocol: the aligned response residual worsened endpoint MAE and did not beat shuffling. |
| Response noise / trajectory length | Five-block response SEM averaged 2.283 kcal mol⁻¹ (median 1.861; maximum 15.087). Noise is material, but longer trajectories were not tested, so a length-based rescue is unproven. |
| λ placement | The same-curve oracle was far better than deployable policies, but held-out-block scoring reversed its advantage and found unstable selected sets; stable placement headroom is not established. |
| Structure prior | The structure-conditioned prior/posterior was not consistently better than generic BQ or direct interpolation, especially for amide, ether, fused-aromatic and alkane sentinels. |
| Cross-fidelity mapping | Not tested. The same-fidelity PIMD2 reconstruction prerequisite failed, so PIMD4/PIMD8 escalation would add degrees of freedom without an established base. |
| Hamiltonian bias | Not a cause of the same-Hamiltonian reconstruction failure. It remains a plausible limit on transfer from short ARROW/PIMD2 responses to experiment. |
| Endpoint labels | The 72-molecule baseline has little headroom, but the aligned-versus-shuffled result is the decisive evidence: no usable molecule-specific endpoint signal was detected. |
| Data availability | No compatible historical dense population existed. The new prospective dense panel resolved that dependency for this bounded test. |
| Implementation failure | Not supported: all 144 new windows passed QC on first attempt, all outputs were finite/complete, and the frozen analysis and unit assertions passed. |

## Compute and failed-work accounting

The campaign generated 144 new windows: 48 calibration and 96 prospective. All passed QC on the first attempt; failed-work cost was zero. New production totaled 720 ps, 288 bead-windows and 900000 nominal bead-steps. Measured simulation wall time was 1.932 h and GPU accounting was 1.932 GPU-h. Arbalest does not expose exact fast/slow force-kernel call counts, so no fabricated force-evaluation number is reported.

## Scientific conclusion

For the inherited short PIMD2 protocol, actual target-molecule response measurements neither improved the experimental endpoint nor enabled the predeclared molecule-conditioned active quadrature to outperform a simple uniform schedule. The result identifies two separable bottlenecks: response residuals were not predictive of experimental error, and the structure-conditioned curve prior/acquisition was not accurate enough to locate the molecule-specific informative windows visible to the oracle.

A future campaign would first require a substantially larger protocol-matched dense-curve training set and independent evidence that its response prior predicts full curves—not only three points—before testing longer trajectories or multi-fidelity escalation. That is a new data-generation program, not a scale-up justified by the present result. Tier-B remains unopened.
