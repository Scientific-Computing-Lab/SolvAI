"""Post-lock replay and size diagnostics for all extension models."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import torch

import tabm_extension as ext
from core import Dataset, SEEDS, dump, outer_partitions, sha, split_inner, checkpoint_predict
from evaluate import load_stress
from refine_representation import RefinedDataset
from deterministic_repair import CorrectedDataset
from fixed_tree_control import tree_arrays


def read(p):
    return json.loads(Path(p).read_text())


def choices(args, d, refined):
    lock = read(args.output / 'extension_locked/external.json')
    yield 'tabm', lock['tabm'], d, 'extension_checkpoints'
    if lock['count_adapted'] is not None:
        yield 'count_adapted', lock['count_adapted'], d, 'extension_checkpoints'
    yield 'embedding_refined', read(args.output / 'representation_locked/external.json'), refined, 'representation_checkpoints'
    yield 'molformer_deterministic', read(args.output / 'deterministic_locked/external.json'), CorrectedDataset(args.release,args.previous), 'deterministic_checkpoints'


def stress(args, d, refined):
    rows = []
    replay_checks = []
    ids = np.arange(56)
    for family, choice, dataset, checkpoint_dir in choices(args, d, refined):
        stressdata, grid = load_stress(args, dataset)
        for response in [True, False]:
            name = family + ('_solvai' if response else '_structure')
            predictions = []
            for seed in SEEDS:
                path = args.output / checkpoint_dir / f'{name}_{seed}.pt'
                cp = torch.load(path, map_location='cpu', weights_only=False)
                assert cp['config'] == choice['candidate'] and cp['response'] == response
                p = ext.replay(stressdata, cp, ids)
                assert np.isfinite(p).all()
                reverse = ext.replay(stressdata, cp, ids[::-1])[::-1]
                delta = float(np.max(np.abs(p - reverse)))
                assert delta < 1e-4
                same = p[(grid.n == 0) & grid.series.isin(['glycine', 'alanine'])]
                assert abs(same[0] - same[1]) < 1e-4
                predictions.append(p)
                replay_checks.append({'model': name, 'seed': seed, 'batch_permutation_max_error': delta,
                                      'checkpoint_sha256': sha(path)})
            x = grid.copy()
            x['model'] = name
            x['prediction'] = np.mean(predictions, axis=0)
            x['seed_std'] = np.std(predictions, axis=0, ddof=1)
            for seed, p in zip(SEEDS, predictions):
                x[f'seed_{seed}'] = p
            rows.append(x)
            print('EXT_STRESS_DONE', name, float(x.prediction.min()), float(x.prediction.max()), flush=True)
    import joblib
    corrected=CorrectedDataset(args.release,args.previous)
    stressdata,grid=load_stress(args,corrected)
    for response in [True,False]:
        name='molformer_fixed_tree_'+('solvai' if response else 'structure');predictions=[]
        x=tree_arrays(stressdata,response)
        for seed in [11,29,47]:
            path=args.output/'deterministic_tree_checkpoints'/f'{name}_{seed}.joblib'
            model=joblib.load(path);p=model.predict(x);predictions.append(p)
            delta=float(np.abs(model.predict(x[::-1])[::-1]-p).max());assert delta<1e-10
            replay_checks.append({'model':name,'seed':seed,'batch_permutation_max_error':delta,'checkpoint_sha256':sha(path)})
        frame=grid.copy();frame['model']=name;frame['prediction']=np.mean(predictions,axis=0)
        frame['seed_std']=np.std(predictions,axis=0,ddof=1)
        for seed,p in zip([11,29,47],predictions):frame[f'seed_{seed}']=p
        rows.append(frame)
    pd.concat(rows, ignore_index=True).to_csv(args.output / 'extension_stress_predictions.csv', index=False)
    dump(args.output / 'extension_stress_replay.json', {'status': 'pass', 'checks': replay_checks,
                                                      'no_experimental_peptide_truth': True})


def size(args, d, refined):
    for family, choice, dataset, _ in choices(args, d, refined):
        outer = np.flatnonzero(dataset.parts != 'external')
        tr, va = split_inner(dataset, outer, 3)
        for response in [True, False]:
            name = family + ('_solvai' if response else '_structure')
            path = args.output / 'extension_size_diagnostic' / f'{name}.json'
            if path.exists():
                continue
            predictions = []
            records = []
            for seed in [11, 29, 47]:
                r, cp, pred = ext.train(dataset, choice['candidate'], tr, np.array([], dtype=int),
                                        seed, response, choice['refit_schedule'])
                predictions.append(pred(va)); records.append(r)
                del cp, pred
                print('EXT_SIZE_FIT', name, seed, flush=True)
            dump(path, {'model': name, 'validation_ids': dataset.frame.id.iloc[va].tolist(),
                        'predictions': np.mean(predictions, axis=0).tolist(),
                        'seed_predictions': [p.tolist() for p in predictions], 'records': records,
                        'selection_used_size_holdout': False,
                        'caveat': 'Exploratory training-only size split; these labels were available during random-fold development.'})


def replay_all(args, d, refined):
    reports = []
    ids = np.flatnonzero(d.parts == 'external')
    specs = [('final_runs', 'checkpoints', d, checkpoint_predict),
             ('extension_final', 'extension_checkpoints', d, ext.replay),
             ('representation_final', 'representation_checkpoints', refined, ext.replay),
             ('deterministic_final','deterministic_checkpoints',CorrectedDataset(args.release,args.previous),ext.replay)]
    for folder, cpdir, dataset, predict in specs:
        for path in sorted((args.output / folder / 'external').glob('*.json')):
            record = read(path)
            cp = torch.load(args.output / cpdir / f'{record["model"]}_{record["seed"]}.pt',
                            map_location='cpu', weights_only=False)
            a = predict(dataset, cp, ids)
            b = np.array(record['predictions'])
            delta = float(np.max(np.abs(a - b)))
            assert delta < 1e-4
            reports.append({'model': record['model'], 'seed': record['seed'], 'max_replay_error': delta})
    expected_families=['resnet','dual','ridge_residual','embedding','graph','tabm','embedding_refined','molformer_deterministic']
    if read(args.output/'extension_locked/external.json')['count_adapted'] is not None:expected_families.append('count_adapted')
    expected={(family+'_'+condition,seed) for family in expected_families for condition in ['solvai','structure'] for seed in SEEDS}
    assert {(r['model'],r['seed']) for r in reports}==expected
    import joblib
    corrected=CorrectedDataset(args.release,args.previous)
    for path in sorted((args.output/'deterministic_tree_final/external').glob('*.json')):
        record=read(path);model=joblib.load(args.output/'deterministic_tree_checkpoints'/f'{record["model"]}_{record["seed"]}.joblib')
        a=model.predict(tree_arrays(corrected,record['with_response'])[ids])
        delta=float(np.abs(a-np.array(record['predictions'])).max());assert delta<1e-10
        reports.append({'model':record['model'],'seed':record['seed'],'max_replay_error':delta})
    expected|={(f'molformer_fixed_tree_{condition}',seed) for condition in ['solvai','structure'] for seed in [11,29,47]}
    assert {(r['model'],r['seed']) for r in reports}==expected
    dump(args.output / 'full_checkpoint_replay.json', {'status': 'pass', 'records': reports})


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['stress', 'size', 'replay'])
    p.add_argument('--release', type=Path, default=Path('RELEASE_ROOT'))
    p.add_argument('--previous', type=Path, default=Path('BASE_INPUTS'))
    p.add_argument('--output', type=Path, default=Path(__file__).parent / 'results')
    args = p.parse_args()
    torch.set_num_threads(2)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    assert (args.output / 'comprehensive_policy.json').exists()
    d, refined = Dataset(args.release, args.previous), RefinedDataset(args.release, args.previous)
    {'stress': stress, 'size': size, 'replay': replay_all}[args.stage](args, d, refined)


if __name__ == '__main__':
    main()
