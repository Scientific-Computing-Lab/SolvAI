"""Check published endpoint results without private contexts or model fitting."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/endpoint_models_20261006'

def read(name):return pd.read_csv(DATA/name,float_precision='round_trip')

def test_endpoint_metrics_reconstruct_from_held_out_rows():
    p=read('combined_predictions.csv');m=read('combined_metrics.csv')
    assert len(p)==11590
    assert not p.duplicated(['model','id']).any()
    for r in m.itertuples():
        g=p.loc[p.model.eq(r.model)]
        g=g.loc[g.part.eq('arrow') if r.cohort=='ARROW-85' else g.part.eq('external')]
        if r.cohort=='Strict-97':g=g.loc[g.strict]
        assert len(g)=={'ARROW-85':85,'External-220':220,'Strict-97':97}[r.cohort]
        error=np.abs(g[r.aggregation]-g.y)
        assert np.isclose(error.mean(),r.mae,atol=1e-12,rtol=0)

def test_paired_contrasts_equal_metric_differences():
    m=read('combined_metrics.csv').set_index(['cohort','model','aggregation'])
    for r in read('combined_comparisons.csv').itertuples():
        a=m.loc[(r.cohort,r.candidate,r.aggregation),'mae']
        b=m.loc[(r.cohort,r.reference,r.aggregation),'mae']
        assert np.isclose(a-b,r.delta,atol=1e-12,rtol=0)
        assert r.low<=r.delta<=r.high

def test_corrected_comparison_preserves_non_molformer_rows():
    old=pd.read_csv(ROOT/'results/journal_20261002/all_predictions.csv',float_precision='round_trip')
    new=read('corrected_original_predictions.csv')
    assert len(new)==2745
    keep=~old.model.str.startswith('molformer_')
    a=old.loc[keep].sort_values(['model','id']).reset_index(drop=True)
    b=new.loc[~new.model.str.startswith('molformer_')].sort_values(['model','id']).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b)
    for cohort,g in new.groupby('part'):
        assert g.id.nunique()=={'arrow':85,'external':220}[cohort]

def test_neutral_framework_and_scientific_figure_formats():
    for name in ['F1_overview','F2_learning_and_evaluation']:
        svg=(ROOT/f'paper/figures/journal/{name}.svg').read_text()
        assert '>AI</text>' in svg
        assert 'data:image/png;base64,' in svg
        assert '/Users/' not in svg and '/home/' not in svg
    for n in range(8,11):
        svg=(ROOT/f'paper/supplementary/figures/Supp_Fig{n}_endpoint_mechanisms.svg').read_text()
        assert 'data:image/png;base64,' in svg
        assert '<text' in svg and '<path' in svg and '#FFFFFF' in svg
        assert '/Users/' not in svg and '/home/' not in svg
    from PIL import Image
    for p in (ROOT/'paper/supplementary/figures/endpoint_components').glob('*.png'):
        a=np.asarray(Image.open(p).convert('RGB'))
        assert (a[:4]==255).all() and (a[-4:]==255).all()
        assert (a[:,:4]==255).all() and (a[:,-4:]==255).all()

def test_public_manifest_and_privacy_boundary():
    manifest=json.loads((DATA/'publication_manifest.json').read_text())
    for name,sha in manifest['files'].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==sha
    for name,record in manifest['implementation'].items():
        p=ROOT/'repro/endpoint_models/implementation'/name
        assert hashlib.sha256(p.read_bytes()).hexdigest()==record['publication_sha256']
    assert not list(DATA.glob('*.joblib')) and not list(DATA.glob('*.npy'))
    for p in [DATA/'selection_specification.json',* (ROOT/'repro/endpoint_models/implementation').glob('*.py')]:
        assert '/home/galoren/' not in p.read_text() and '/Users/galoren/' not in p.read_text()
