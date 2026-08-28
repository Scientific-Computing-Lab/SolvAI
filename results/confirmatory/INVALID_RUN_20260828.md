# Invalid environment diagnostic

The initial diagnostic call to `scripts/run_confirmatory_endpoint.py --mode primary`
used `/home/galoren/anaconda3/bin/python` with scikit-learn 1.6.1. The frozen release
environment requires scikit-learn 1.7.2 and is invoked through `uv run` from the
parent workspace. The preliminary prediction, metric, comparison and metadata files
were moved into this directory under the `invalid_sklearn_1_6_1_` prefix and are not
used in any confirmatory result or conclusion.
