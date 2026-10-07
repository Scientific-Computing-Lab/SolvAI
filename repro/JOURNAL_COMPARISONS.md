# Reproducing the journal molecular-model comparisons

The October 2026 comparison retains the original ARROW-85 fixed folds and the
external 220-molecule cohort (including the fixed strict subset of 97). It is a
retrospective extension, not a new prospective test. All families and seeds are
reported. The original SolvAI architecture and response surrogates are unchanged.

## Inputs and environments

The original data workspace must contain the processed source tables used by
`scripts/confirmatory_common.py`. The release must contain the standardized-exclusion
surrogate predictions and the original external-cohort records. Follow
`repro/DATA_PROVENANCE.md` for source acquisition and licensing. Raw third-party
training labels and generated feature caches remain local; the protocol, input hashes,
test predictions, validation traces and aggregate results are distributable outputs.
Published input-hash keys are repository-relative paths; the digest values are
unchanged from the original execution manifest.

Graph fits use Chemprop 2.2.4, PyTorch 2.5.1 and scikit-learn 1.7.2. Uni-Mol embedding
generation uses the original separate Uni-Mol environment. MoLFormer uses the pinned
revision in `run_journal_baselines.py`, loaded from the existing verified model cache.
The published CheMeleon weights must be available as `.chemprop/chemeleon_mp.pt`
under the executing user's home. No API key or paid model inference is required.

## Sequence

The commands below document the original fixed-configuration study. Its MoLFormer
cache was subsequently found to retain stochastic attention projections. The
current paper uses the deterministic correction in `repro/endpoint_models/` and
`results/endpoint_models_20261006/corrected_original_*` instead of those legacy
MoLFormer rows. Do not regenerate a current claim from the legacy cache. For the
published tables and figures, use `python scripts/make_endpoint_figures.py` followed
by `python scripts/make_local_journal_figures.py`; neither command refits models.

Set `SOLVAI_WORKSPACE`, `SOLVAI_RELEASE` and `SOLVAI_RESULTS` to the original data
workspace, this repository and a results directory. Use the appropriate environment's
Python executable for each command:

```sh
python scripts/run_journal_baselines.py prepare --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/run_journal_baselines.py embeddings --encoder molformer --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/run_journal_baselines.py embeddings --encoder unimol --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/run_journal_baselines.py computed --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/run_journal_baselines.py trees --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/run_journal_baselines.py graphs --workspace "$SOLVAI_WORKSPACE" --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
python scripts/summarize_journal_baselines.py --release "$SOLVAI_RELEASE" --output "$SOLVAI_RESULTS"
make journal-paper
```

The two graph families may run independently with `--graph-variants dmpnn` and
`--graph-variants chemeleon`; do not launch duplicate fits for the same family into
one output directory. Completed seed runs are reused. Graph caching changes only
feature-construction cost, not the graph representation.

## Verification and interpretation

The summary requires all eight primary model conditions and the computation-only
control, all five ARROW folds and all 220 external predictions.
It checks seed averaging, molecule identities and exact target alignment. Refit
structural and SolvAI controls must reproduce the released external predictions
to a maximum absolute difference below `1e-10` kcal/mol before any new comparison
is accepted. CSV training labels must use round-trip floating-point parsing; tiny
target perturbations can change tree split tie-breaking.

The graph-model epoch count is selected on an inner validation split, followed by
a fresh fit using all outer-training labels. Test labels never select epochs or
hyperparameters. Frozen representation comparisons use the same ExtraTrees learner
and compare absence/presence of the 15 response descriptors. The computation-only
control is a separately specified source-attribution analysis, first on ARROW and
then extended without retuning to the external cohort. The external OpenFF/GBn2
components must reconstruct the original corrected response sums before fitting.

The external error strata are descriptive; they are not a calibrated applicability
domain. Molecular membership in foundation-model pretraining is not audited. No
claim of an exhaustive architecture or hyperparameter ranking follows from these
fixed configurations. See `repro/JOURNAL_FIGURES.md` for the current eight main and
twelve supplementary figures, including the later endpoint study. The main text, references and Supplementary
Information are compiled together in `paper/main.pdf` by `make journal-paper`.
