# Freecurve repository recovery

## Outcome

Authenticated GitHub access recovered **all 40 repositories** in the supplied
organization manifest into `repositories_recovered/`. This includes the private
repositories that were absent from the original six-repository package:
`ARBALEST_EXAMPLES`, `FF-NN`, `nndb-extracted`,
`normalizing_flows-extracted`, `arbalest-app-data`, `FF-Tools`, and the remaining
organization support repositories.

The recovery used a temporary, isolated GitHub configuration. The pre-existing
global GitHub identity was not changed, and no credential is stored in this
workspace, source code, Git remote URL, report, or catalog. Recovered remotes use
normal `https://github.com/freecurvelabs/...` URLs.

After the history audit, the temporary authentication directory and all
credential-bearing Git worker processes were removed. The access token pasted
into the conversation must nevertheless be revoked/rotated because chat text is
not a safe credential channel.

## Scientifically relevant additions

- The current `arbalest` default branch contains active 15-window TI/PIMD
  campaigns for water, methane, propane, pyrrole, butane, and cyclohexane,
  including committed comparison plots and precise run manifests. The raw run
  directories referenced by the manifests are not committed; recovery through
  GitHub artifacts/history is tracked separately from image digitization.
- `FF-NN` and `normalizing_flows-extracted` contain large component-resolved
  water-dimer training collections and ARROW/GAFF paired interaction energies.
- `FF-Tools` contains 22K acetamide-water structures, directly relevant to a
  current high-error chemistry family.
- `arbalest-app-data` contains NMA-water, butane-water, propane, and other
  trajectory/structure assets that were not present in the initial package.
- `scripts/arbalest/ene-temp` contains a complete 21-window, component-resolved
  TI campaign. Its molecular identity is not established yet, so it is excluded
  from structure-conditioned training until provenance is resolved.
- `Sampling` contains 400 water-only HREX snapshots. These are useful background
  solvent data but are not molecule-aligned hydration supervision.

The follow-up all-history object/path scan found **3,358 physics-related Git
paths** across the recovered organization. Historical `FF-NN` objects add
NMA-water, NMA-NMA, and water-dimer `energyStat` archives; `Simulations` adds one
complete toluene TI regression fixture. No hidden second production-scale,
molecule-resolved hydration campaign with reusable lambda observables was found.
The recovered intermolecular archives are useful for representation and
interaction-energy supervision, but they cannot supply the missing
per-molecule ARROW/PIMD response labels.

## Reproducibility

- Source manifest: `REPOSITORY_MANIFEST.csv`
- Recovery catalog: `data/catalog/freecurve_repository_recovery.parquet`
- Human-readable catalog: `data/catalog/freecurve_repository_recovery.csv`
- Recovery root: `repositories_recovered/`
- Recovered worktrees verified clean: **40/40**
- New production PIMD8 labels generated during recovery: **0**

Only labels with a resolved molecule, method, units, sign convention, and source
commit will enter the physics-teacher table.
