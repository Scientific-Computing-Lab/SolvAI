# Active SolvAI final handoff

> **Subsequent qualification:** the original same-curve oracle values below do
> not independently establish stable placement headroom. Held-out-block scoring
> reversed their advantage, and the independent-replica resolution campaign was
> blocked by its prospective power gate before simulation. The immutable no-go
> decision itself is unchanged.

## Repository state

- Parent SolvAI reference: `531f6cfd21e319c951b461c9ef24fa754790f91d`
- Isolated branch: `active-solvai`
- Blueprint SHA-256: `957bd6b3a8780a921244069c73a164cf1bce49aefd51fd6bc88e3144f7c86372`
- Program decision: **NO-GO**
- Tier-B status: unopened; no blind labels or predictions were accessed

## Decisive results

1. **Experimental endpoint:** actual-minus-predicted PIMD2 residuals changed
   five-partition MAE from 0.186374 to 0.189856 kcal mol⁻¹. The paired change
   was +0.003482 (95% CI +0.000960, +0.006093), and aligned residuals did not
   beat shuffled controls.
2. **Prospective dense reconstruction:** active structure-conditioned Bayesian
   quadrature had MAEs of 1.701 and 1.608 kcal mol⁻¹ at five and seven windows;
   the strongest comparator, uniform direct integration, had 1.153 and 1.092.
3. **Placement diagnostic:** the non-deployable same-curve oracle reached 0.337
   and 0.068 kcal mol⁻¹ at the same budgets. Subsequent held-out-block scoring
   reversed this apparent advantage and found unstable selected sets, so stable
   molecule-specific placement headroom is not established.
4. **Escalation:** Direction C was prospectively conditional on Direction B and
   was not launched. No PIMD4/PIMD8 scale-up, new Hamiltonian, learned potential
   or Tier-B campaign was opened.

## Compute accounting

The dense campaign generated 144 new 5-ps PIMD2 windows: 48 calibration and 96
prospective. All passed quality control on the first attempt. Total production
was 720 ps, 288 bead-windows, 900,000 nominal bead-steps, 6,954.3 seconds of
measured simulation wall time and 1.932 GPU-h. Failed-work cost was zero. Exact
fast/slow force-kernel evaluations are unavailable from Arbalest and were not
invented.

## Canonical deliverables

- Scientific report: `reports/ACTIVE_SOLVAI_FINAL_REPORT.md`
- Phase 1 report: `reports/PHASE1_ACTUAL_OBSERVATION_GATE.md`
- Prospective replay report: `reports/DENSE_SENTINEL_REPLAY.md`
- Canonical metrics: `results/phase2/active_solvai_final_metrics.json`
- Molecule diagnostics: `results/phase2/dense_failure_diagnostics.csv`
- Molecule-level replay: `results/phase2/dense_replay_predictions.parquet`
- Main manuscript: `paper/main.pdf`
- Supplement: `paper/supplementary/supplementary.pdf`
- Immutable ledger: `runs/ledger.jsonl`
- File manifest: `MANIFEST.json`

## Scale decision

No positive scale request is justified. The most direct future requirement is
a substantially larger, protocol-matched dense-curve corpus that improves
structure-to-curve prediction and prospectively identifies informative λ
locations. This would be a new program requiring its own freeze and budget; it
is not a continuation authorized by the present negative result.
