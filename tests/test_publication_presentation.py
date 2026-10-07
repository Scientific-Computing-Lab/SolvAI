"""Publication labels must match the scientific definitions and saved values."""
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def test_manuscript_definitions_and_signed_reduction():
    main = (ROOT/'paper/main.tex').read_text()
    supp = (ROOT/'paper/supplementary/supplementary.tex').read_text()
    assert r'\linenumbers' not in main
    assert '1,280 external' not in main+supp
    assert 'maximum allowed train--test' not in main
    assert 'eight-bead PIMD (PIMD8)' in main
    assert 'four-bead PIMD (PIMD4)' in main
    assert 'nuclear quantum effect (NQE)' in main
    assert 'Root mean square error (RMSE)' in main
    assert r'reduces MAE by \SolvAIReduction{}' in main
    macros = (ROOT/'paper/tables/metrics_macros.tex').read_text()
    metrics = json.loads((ROOT/'results/paper_metrics.json').read_text())
    delta = next(r['difference'] for r in metrics['paired_confirmatory']
                 if r['analysis']=='primary_F_full_solvai')
    value = float(re.search(r'\\newcommand\{\\SolvAIReduction\}\{([^}]+)\}',macros)[1])
    assert value == round(abs(delta),3)
    assert r'\ensuremath{-0.101}' in macros


def test_supplementary_data_file_map_and_units():
    supp = (ROOT/'paper/supplementary/supplementary.tex').read_text()
    for name in re.findall(r'\\nolinkurl\{(Supplementary_Data_[^}]+)\}',supp):
        assert (ROOT/'paper/supplementary_data'/name).is_file(), name
    for i in range(1,6):
        assert f'Supplementary Data {i}:' in supp
    from scripts import make_figures
    from unittest.mock import patch
    import pandas as pd
    data=pd.read_parquet(ROOT/'results/confirmatory/standardized_exclusion_endpoint_predictions.parquet')
    primary=data.loc[data.partition.eq('standardized_exclusion_primary')]
    separation=pd.read_csv(ROOT/'results/confirmatory/standardized_exclusion_global_separation_metrics.csv')
    with patch.object(make_figures,'save') as saved:
        make_figures.supp_fig5_extrapolation(primary,separation)
        figure=saved.call_args.args[0]
        assert 'mol$^{-1}$' in figure.axes[0].get_xlabel()
        assert 'NN < 0.70' in [t.get_text() for t in figure.axes[1].texts]
        make_figures.plt.close(figure)


def test_lambda_redesign_keeps_every_frozen_value():
    from scripts import make_figures
    from unittest.mock import patch
    import numpy as np
    metrics=json.loads((ROOT/'results/paper_metrics.json').read_text())
    expected=list(metrics['multilambda']['method_mae_kcal_mol'].values())
    with patch.object(make_figures,'save') as saved:
        make_figures.supp_fig4_lambda(metrics)
        fig=saved.call_args.args[0]
        assert len(fig.axes)==3
        np.testing.assert_allclose([b.get_width() for b in fig.axes[1].patches],expected[:4],rtol=0,atol=0)
        assert float(fig.axes[2].collections[0].get_offsets()[0,1])==expected[4]
        assert len(fig.axes[0].lines)==3
        make_figures.plt.close(fig)
