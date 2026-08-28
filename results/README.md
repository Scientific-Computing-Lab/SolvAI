# Results provenance

`paper_metrics.json`, `paper_metrics.csv`, `final_metrics.json` and
`model_comparison.csv` are synchronized canonical manuscript outputs. They are
regenerated from molecule-level files in `confirmatory/`.

`oof_predictions.parquet` and the `predictions/`, `robustness/` and `ablations/`
directories preserve the exploratory campaign. They include the historical 0.19705
exact-connectivity analysis and must not be used as the publication headline. Those
values remain explicitly labelled under `historical_campaign` in the canonical JSON.
