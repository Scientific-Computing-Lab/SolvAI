# Endpoint architectures and reproducibility

The manuscript's extended endpoint study is retrospective. The original
ExtraTrees release is unchanged. No new prospective validation or globally
optimal architecture is claimed.

## Frozen-output reproduction

`results/endpoint_models_20261006/` contains molecule-level test and unlabelled-series
predictions, all family/ablation metrics, common-three-seed comparisons, paired
intervals, size-diagnostic aggregates, native interval coverage and configuration
choices. `scripts/make_endpoint_figures.py` rebuilds the numeric tables and native
vector figures locally. It also replaces the earlier MoLFormer comparison with
deterministically encoded predictions without altering the old source archive.

`combined_*` retains legacy stochastic-cache variants for an audit trail; rows with
`legacy_stochastic_molformer` or nonzero `legacy_stochastic_molformer_rows` are
NOT valid performance claims. Published figures filter them out. The corrected
MoLFormer is identified by `molformer_deterministic` or `molformer_fixed_tree`.
`capacity_adapted`, `count_adapted` and `selected_*` are selection policies, not
new network architectures. Strict-97 is nested within External-220.

## Fitting and restricted dependencies

`implementation/` archives the actual architecture, preprocessing, fitting and
evaluation code. Only machine-specific path strings and introductory descriptions
were relocated; the manifest records both original and exported hashes. This is
implementation source, not a bundled training dataset or one-command deployment.
Machine path placeholders `RELEASE_ROOT`, `BASE_INPUTS` and `STUDY_ROOT` must point
to acquired inputs and an isolated working directory before rerunning.

Required inputs are the original released surrogates, endpoint identities/features
and source labels acquired under their original terms, corrected MoLFormer
encodings and the recorded splits/configurations. See `repro/DATA_PROVENANCE.md`.
The immutable experiment-input hashes and candidate/selection specifications are
in `selection_specification.json`. Training pools and fitted TabPFN contexts are
not redistributed: the latter embed training rows and labels. In particular,
the size test reuses 205 labels from the development pool; its aggregate metrics
are released, but restricted source-label rows must be obtained from their sources.

Runtime used Python 3.11.15, PyTorch 2.5.1+cu121 and Chemprop 2.2.4, TabM 0.0.3,
rtdl-num-embeddings 0.0.12, and isolated TabPFN 9.1.0. Third-party packages are
installed from their original distributions, not vendored here. TabPFN-3.5 full
and Fast checkpoints require the Prior Labs non-commercial license and revision
06bf2ba35c80a92a3b9abb436b99cf49e7a0365e. No foundation-model fine-tuning was done.
The full experiment's hashed audit records remain retained by the authors.

## Interpretation

Configurations, schedules and feature recipes were chosen inside each outer
training pool; matched no-response controls inherit these recipes. Existing
cohorts had already been examined in earlier development, so locked final choices
do not make this a new prospective test. Bootstrap intervals use 100,000 paired
molecule resamples and are conditional/descriptive, without multiplicity adjustment.
Homologous peptide predictions have no reference values; smoothness, magnitude
and absence of saturation are not measures of physical accuracy.
