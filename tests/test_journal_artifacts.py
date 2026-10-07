"""Lightweight checks of the frozen journal outputs; no model fitting required."""
import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/journal_20261002"
MODELS = {
    "structure", "solvai", "computed_only", "dmpnn", "chemeleon",
    "molformer_structure", "molformer_solvai", "unimol_structure", "unimol_solvai",
}


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_journal_predictions_and_metrics_match():
    predictions = read_csv(RESULTS / "all_predictions.csv")
    assert len(predictions) == 2745
    assert {row["model"] for row in predictions} == MODELS
    assert len({(row["model"], row["id"]) for row in predictions}) == len(predictions)
    for row in predictions:
        mean = sum(float(row[f"seed_{seed}"]) for seed in (11, 29, 47)) / 3
        assert abs(mean - float(row["prediction"])) < 1e-12
    for metric in read_csv(RESULTS / "metrics.csv"):
        rows = [row for row in predictions if row["model"] == metric["model"]]
        if metric["cohort"] == "ARROW-85":
            rows = [row for row in rows if row["part"] == "arrow"]
            assert len(rows) == 85
        else:
            rows = [row for row in rows if row["part"] == "external"]
            if metric["cohort"] == "Strict-97":
                rows = [row for row in rows if row["strict"].lower() == "true"]
                assert len(rows) == 97
            else:
                assert len(rows) == 220
        mae = sum(abs(float(row["y"]) - float(row["prediction"])) for row in rows) / len(rows)
        assert abs(mae - float(metric["mae"])) < 1e-12


def test_journal_audit_and_script_identity():
    audit = json.loads((RESULTS / "audit.json").read_text())
    validation = json.loads((RESULTS / "validation.json").read_text())
    assert audit["status"] == "pass" and audit["graph_seed_runs"] == 36
    assert audit["outer_test_labels_entered_fitting_or_epoch_selection"] is False
    for report, field, script in [
        (audit, "audit_script_sha256", "audit_journal_baselines.py"),
        (validation, "analysis_script_sha256", "summarize_journal_baselines.py"),
    ]:
        assert report[field] == hashlib.sha256((ROOT / "scripts" / script).read_bytes()).hexdigest()
    assert all(value < 1e-10 for value in validation["control_replay_max_absolute_difference"].values())
    assert len(list((RESULTS / "predictions").glob("*.csv"))) == 54


def test_journal_figures_use_vector_assets():
    main = (ROOT / "paper/main.tex").read_text()
    supplementary = (ROOT / "paper/supplementary/supplementary.tex").read_text()
    main += (ROOT / "paper/endpoint_results.tex").read_text()
    main += (ROOT / "paper/endpoint_size_results.tex").read_text()
    assert ".png}" not in main + supplementary
    assert ".jpg}" not in main + supplementary
    names=['F1_overview','F2_learning_and_evaluation','F3_endpoint_comparison',
           'F4_matched_evidence','F5_modern_comparisons','F6_generalization',
           'F7_source_composition','F8_size_and_stress']
    folder=ROOT/'paper/figures/journal'
    manifest=json.loads((folder/'sources.json').read_text())
    assert manifest['main_figures']==names
    assert manifest['source_descriptor_counts']==[1,5,1,1,1,6]
    assert manifest['computed_only_descriptor_count']==10
    for name in names:
        assert f'{{{name}.pdf}}' in main
        assert (folder/f'{name}.pdf').is_file()
        svg=(folder/f'{name}.svg').read_text()
        assert bool(re.search(r'<(?:\w+:)?image\b',svg)) == (name == 'F1_overview')
        assert '<text' in svg
    for field in ['source_sha256','builder_sha256']:
        for path,digest in manifest[field].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    for number,stem in [(6,'evaluation_protocol'),(7,'conformer_targets')]:
        name=f'Supp_Fig{number}_{stem}'
        assert f'{{{name}.pdf}}' in supplementary
        assert f'Supplementary Fig.~{number}' in main
        svg=(ROOT/f'paper/supplementary/figures/{name}.svg').read_text()
        assert not re.search(r'<(?:\w+:)?image\b',svg)


def test_supplementary_screen_baselines_are_campaign_specific():
    from unittest.mock import patch
    from scripts import make_figures

    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    campaigns = make_figures.supp_fig3_campaigns(metrics)
    assert [len(frame) for _, _, frame in campaigns] == [5, 2]
    source_base = float(read_csv(ROOT / "results/ablations/smd_plus_confsolv_response_screen.csv")[0]["mae"])
    lambda_base = float(read_csv(ROOT / "results/ablations/multilambda_metrics.csv")[0]["mae"])
    assert abs(campaigns[0][1] - source_base) < 1e-12
    assert abs(campaigns[1][1] - lambda_base) < 1e-12
    assert source_base != lambda_base
    expected_source_files = {
        "OpenFE diagnostics": "openfe_screen.csv",
        "MLFF hierarchy": "mlff_screen.csv",
        "DES370K response": "des370k_screen.csv",
        "ConfSolv graph latent": "smd_confsolv_graph_embedding_screen.csv",
        "ConfSolv FFN latent": "smd_confsolv_ffn_embedding_screen.csv",
    }
    for name, source in expected_source_files.items():
        assert abs(campaigns[0][2][name] - float(read_csv(ROOT / "results/ablations" / source)[0]["mae"])) < 1e-10
    assert set(campaigns[1][2].index) == {"PIMD2 lambda response", "Classical/NQE/PIMD"}
    with patch.object(make_figures, "save") as save:
        make_figures.supp_fig3_alternatives(metrics)
        figure = save.call_args.args[0]
        assert len(figure.axes) == 2
        for axis, (_, baseline, frame) in zip(figure.axes, campaigns, strict=True):
            assert len(axis.lines) == 1
            assert list(axis.lines[0].get_xdata()) == [baseline, baseline]
            assert list(axis.collections[1].get_offsets()[:, 0]) == list(frame.values)
            assert [label.get_text() for label in axis.get_yticklabels()] == list(frame.index)
        make_figures.plt.close(figure)


def test_supplementary_residual_panels_do_not_overlap_or_change_data():
    from unittest.mock import patch
    import numpy as np
    import pandas as pd
    from scripts import make_figures

    data = pd.read_parquet(ROOT / 'results/confirmatory/standardized_exclusion_endpoint_predictions.parquet')
    primary = data.loc[data.partition.eq('standardized_exclusion_primary')]
    full = primary.loc[primary.method.eq('F_full_solvai')].sort_values('molecule_id')
    base = primary.loc[primary.method.eq('A_structure_only')].sort_values('molecule_id')
    with patch.object(make_figures, 'save') as save:
        make_figures.supp_fig1_residuals(primary)
        figure = save.call_args.args[0]
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        boxes = [axis.get_tightbbox(renderer) for axis in figure.axes]
        assert all(a.x1+5 < b.x0 for a, b in zip(boxes, boxes[1:]))
        bounds = [axis.get_window_extent(renderer) for axis in figure.axes]
        assert max(b.y0 for b in bounds)-min(b.y0 for b in bounds) < .5
        assert max(b.y1 for b in bounds)-min(b.y1 for b in bounds) < .5
        np.testing.assert_array_equal(figure.axes[0].collections[0].get_offsets(),
                                      np.column_stack([full.y_true, full.y_pred]))
        np.testing.assert_array_equal(figure.axes[1].collections[0].get_offsets(),
                                      np.column_stack([full.y_true, full.y_pred-full.y_true]))
        delta = full.absolute_error.to_numpy()-base.absolute_error.to_numpy()
        counts, _ = np.histogram(delta, bins=np.linspace(delta.min(), delta.max(), 17))
        np.testing.assert_array_equal([bar.get_height() for bar in figure.axes[2].patches], counts)
        assert counts.sum() == 85
        make_figures.plt.close(figure)


def test_journal_tree_edges_meet_nodes_and_repeated_scores_are_preserved():
    from unittest.mock import patch
    import sys
    import numpy as np
    from matplotlib.collections import PathCollection
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    sys.path.insert(0, str(ROOT / 'scripts'))
    import make_local_journal_figures as journal
    art = journal.Art(100, 'Geometry test')
    edges, nodes = [], []
    with patch.object(art, 'line', side_effect=lambda x1,y1,x2,y2,*a,**kw: edges.append(((x1,y1),(x2,y2)))), \
         patch.object(art, 'circle', side_effect=lambda x,y,*a,**kw: nodes.append((x,y))):
        art.trees(0,0,120,52)
    assert len(edges) == 18 and len(nodes) == len(set(nodes)) == 21
    assert all(start in nodes and end in nodes for start,end in edges)

    plots = {}
    def collect(self, figure, x, y, w, h, prefix):
        plots[prefix] = figure
    data = journal.load_data()
    with patch.object(journal.Art, 'embed_plot', collect), patch.object(journal.Art, 'save'):
        journal.matched(data, None)
    repeat = data['repeats'].groupby(['repeat','method']).absolute_error.mean().unstack()
    axes = plots['repeats'].axes[0]
    scatters = [c for c in axes.collections if isinstance(c, PathCollection)]
    assert len(scatters) == 2
    for scatter, method in zip(scatters, ['A_structure_only','F_full_solvai']):
        np.testing.assert_array_equal(scatter.get_offsets()[:,1], repeat[method])
        assert len(set(scatter.get_offsets()[:,0])) == 5
    for key in ['progression','repeats']:
        figure = plots[key]
        FigureCanvasAgg(figure)
        figure.canvas.draw()
        axis = figure.axes[0]
        label = next(text for text in axis.texts if text.get_text() == 'PIMD8')
        box = label.get_window_extent(figure.canvas.get_renderer()).padded(3)
        assert box.y0 > axis.bbox.y0
        for scatter in [c for c in axis.collections if isinstance(c, PathCollection)]:
            positions = axis.transData.transform(scatter.get_offsets())
            assert not any(box.contains(x,y) for x,y in positions)


def test_native_tree_port_geometry_and_endpoint_mechanisms():
    from unittest.mock import patch
    import sys
    sys.path.insert(0,str(ROOT/'scripts'))
    import journal_diagram as d
    c=d.Art(120,'Native tree geometry')
    edges=[];circles=[];leaves=[]
    with patch.object(c,'line',side_effect=lambda x1,y1,x2,y2,*a,**k:edges.append(((x1,y1),(x2,y2)))), \
         patch.object(c,'circle',side_effect=lambda x,y,*a,**k:circles.append((x,y))), \
         patch.object(c,'rect',side_effect=lambda x,y,w,h,*a,**k:leaves.append((x+w/2,y+h/2))):
        root,output=d.tree(c,5,10,100,70,selected=2)
    nodes=circles+leaves
    assert len(edges)==6 and len(nodes)==len(set(nodes))==7
    assert all(a in nodes and b in nodes for a,b in edges)
    assert root==circles[0] and output==(leaves[2][0],leaves[2][1]+3.5)
    for n in [8,9,10]:
        svg=(ROOT/f'paper/supplementary/figures/Supp_Fig{n}_endpoint_mechanisms.svg').read_text()
        assert '<image' not in svg and 'data:image/' not in svg
    svg=(ROOT/'paper/supplementary/figures/Supp_Fig10_endpoint_mechanisms.svg').read_text()
    assert 'diag(sᵢ) W diag(rᵢ) x' in svg
