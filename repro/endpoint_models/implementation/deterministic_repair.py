"""Reproducible MoLFormer representation repair for the endpoint comparison."""
from pathlib import Path
import argparse
import json
import numpy as np
import torch
import tabm_extension as ext
import refine_representation as rep
from core import SEEDS, candidates as original_candidates, dump, outer_partitions, schedule_from_records, sha, split_inner

HERE = Path(__file__).parent


def read(p): return json.loads(Path(p).read_text())


def hashes():
    result = {**rep.hashes()}
    for name in ['deterministic_repair.py', 'encode_deterministic_molformer.py', 'DETERMINISTIC_REPAIR.md',
                 'results/deterministic_encoding.json', 'results/molformer_deterministic.npy',
                 'results/stress_molformer_deterministic.npy']:
        result[name] = sha(HERE / name)
    return result


class CorrectedDataset(rep.RefinedDataset):
    def __init__(self, release, previous):
        super().__init__(release, previous)
        freeze = read(HERE / 'results/deterministic_encoding.json')
        assert freeze['status'] == 'pass' and freeze['deterministic_eval']
        for name, digest in freeze['arrays'].items(): assert sha(HERE / 'results' / name) == digest
        self.embeddings['molformer'] = np.load(HERE / 'results/molformer_deterministic.npy')
        assert self.embeddings['molformer'].shape == (len(self.frame), 768)


def candidates():
    result = []
    for source in [c for c in original_candidates() if c.get('encoder') == 'molformer']:
        for recipe in ['binary', 'raw_count', 'log_count', 'no_fingerprint']:
            result.append({**source, 'id': f'd{len(result):02d}',
                           'representation': 'binary_log' if recipe == 'binary' else 'count_log',
                           'count_scaling': recipe, 'input_recipe': recipe,
                           'encoder_pipeline': 'deterministic', 'source_config_id': source['id']})
    return result


def inner(args, d, part, outer, c, split, seed, tag='ordinary'):
    path = args.output / 'deterministic_inner' / part / f'{c["id"]}_{tag}_split{split}_seed{seed}.json'
    if path.exists():
        r = read(path); assert r['config'] == c and r['deterministic_hashes'] == hashes(); return r
    print('DET_START', part, c['id'], c['input_recipe'], tag, split, seed, flush=True)
    tr, va = split_inner(d, outer, split)
    r, cp, pred = ext.train(d, c, tr, va, seed)
    r.update(partition=part, split=split, tag=tag, deterministic_hashes=hashes())
    dump(path, r); del cp, pred
    print('DET_DONE', part, c['id'], tag, split, seed, round(r['best_score'], 5), r['best_epoch'], flush=True)
    return r


def search(args, d):
    assert read(args.output / 'deterministic_protocol.json')['hashes'] == hashes()
    for part, outer, test in outer_partitions(d):
        if args.parts and part not in args.parts: continue
        if (args.output / 'deterministic_locked' / f'{part}.json').exists(): continue
        screen = [(inner(args, d, part, outer, c, 0, 11)['best_score'], c) for c in candidates()]
        repeated = []
        for _, c in sorted(screen, key=lambda r: (r[0], r[1]['id']))[:2]:
            rows = [inner(args, d, part, outer, c, j, s) for j in range(3) for s in [11, 29]]
            repeated.append((float(np.mean([r['best_score'] for r in rows])), c))
        _, c = min(repeated, key=lambda r: (r[0], r[1]['id']))
        longer = {**c, 'cap': 1200, 'schedule': 'cosine'}
        normal = [inner(args, d, part, outer, c, j, 47) for j in range(3)]
        long = [inner(args, d, part, outer, longer, j, 47, 'long') for j in range(3)]
        diffs = np.array([b['best_score'] - a['best_score'] for a, b in zip(normal, long)])
        use_long = bool(diffs.mean() < -.005 and (diffs < 0).sum() >= 2)
        if use_long:
            c = longer
            rows = long + [inner(args, d, part, outer, c, j, s, 'long') for j in range(3) for s in [11, 29]]
        else:
            rows = [inner(args, d, part, outer, c, j, s) for j in range(3) for s in [11, 29, 47]]
        dump(args.output / 'deterministic_locked' / f'{part}.json', {
            'candidate': c, 'mean_score': float(np.mean([r['best_score'] for r in rows])),
            'use_long': use_long, 'long_minus_normal': diffs.tolist(),
            'refit_schedule': schedule_from_records(rows), 'chosen_epochs': [r['best_epoch'] for r in rows],
            'deterministic_hashes': hashes(), 'test_labels_used': False})
        print('DET_LOCKED', part, c['id'], flush=True)


def final(args, d):
    assert (args.output / 'comprehensive_policy.json').exists()
    for part, _, _ in outer_partitions(d): assert (args.output / 'deterministic_locked' / f'{part}.json').exists()
    for part, outer, test in outer_partitions(d):
        if args.parts and part not in args.parts: continue
        choice = read(args.output / 'deterministic_locked' / f'{part}.json')
        for response in [True, False]:
            name = 'molformer_deterministic_' + ('solvai' if response else 'structure'); values = []
            for seed in SEEDS:
                path = args.output / 'deterministic_final' / part / f'{name}_{seed}.json'
                if path.exists(): values.append(np.array(read(path)['predictions'])); continue
                r, cp, pred = ext.train(d, choice['candidate'], outer, np.array([], dtype=int), seed,
                                        response, choice['refit_schedule'], part == 'external')
                y = pred(test)
                r.update(model=name, partition=part, test_ids=d.frame.id.iloc[test].tolist(),
                         predictions=y.tolist(), deterministic_hashes=hashes())
                if cp is not None:
                    cp_path = args.output / 'deterministic_checkpoints' / f'{name}_{seed}.pt'
                    cp_path.parent.mkdir(exist_ok=True); torch.save(cp, cp_path)
                    delta = float(np.max(np.abs(ext.replay(d, cp, test) - y))); assert delta < 1e-4
                    r.update(checkpoint_sha256=sha(cp_path), replay_max_error=delta)
                dump(path, r); values.append(y); del cp, pred
                print('DET_REFIT', part, name, seed, flush=True)
            x = d.frame.iloc[test].copy(); x['model'] = name; x['prediction'] = np.mean(values, axis=0)
            for seed, y in zip(SEEDS, values): x[f'seed_{seed}'] = y
            path = args.output / 'deterministic_predictions' / f'{name}_{part}.csv'
            path.parent.mkdir(exist_ok=True); x.to_csv(path, index=False)


def tree_arrays(d, response):
    blocks = [d.struct, d.embeddings['molformer']]
    if response: blocks.append(d.response)
    return np.column_stack(blocks)


def fixed_tree(args, d):
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline
    import joblib
    assert (args.output / 'comprehensive_policy.json').exists()
    for part, outer, test in outer_partitions(d):
        if args.parts and part not in args.parts: continue
        for response in [True, False]:
            name = 'molformer_fixed_tree_' + ('solvai' if response else 'structure'); values = []
            x = tree_arrays(d, response)
            for seed in [11, 29, 47]:
                path = args.output / 'deterministic_tree_final' / part / f'{name}_{seed}.json'
                if path.exists(): values.append(np.array(read(path)['predictions'])); continue
                model = make_pipeline(SimpleImputer(strategy='median', add_indicator=True),
                    ExtraTreesRegressor(n_estimators=360, min_samples_leaf=2, max_features=.7, random_state=seed, n_jobs=2))
                model.fit(x[outer], d.y[outer], extratreesregressor__sample_weight=d.weights[outer])
                y = model.predict(x[test]); values.append(y)
                r = {'model': name, 'seed': seed, 'partition': part, 'with_response': response,
                     'train_ids': d.frame.id.iloc[outer].tolist(), 'validation_ids': [],
                     'test_ids': d.frame.id.iloc[test].tolist(), 'predictions': y.tolist(),
                     'fixed_config': {'n_estimators': 360, 'min_samples_leaf': 2, 'max_features': .7},
                     'deterministic_hashes': hashes()}
                if part == 'external':
                    cp_path = args.output / 'deterministic_tree_checkpoints' / f'{name}_{seed}.joblib'
                    cp_path.parent.mkdir(exist_ok=True); joblib.dump(model, cp_path, compress=3)
                    delta = float(np.abs(joblib.load(cp_path).predict(x[test]) - y).max()); assert delta < 1e-10
                    r.update(checkpoint_sha256=sha(cp_path), replay_max_error=delta)
                dump(path, r)
                print('DET_TREE_REFIT', part, name, seed, flush=True)
            frame = d.frame.iloc[test].copy(); frame['model'] = name; frame['prediction'] = np.mean(values, axis=0)
            for seed, y in zip([11, 29, 47], values): frame[f'seed_{seed}'] = y
            path = args.output / 'deterministic_tree_predictions' / f'{name}_{part}.csv'
            path.parent.mkdir(exist_ok=True); frame.to_csv(path, index=False)


def prepare(args, d):
    for folder in ['predictions', 'extension_predictions', 'representation_predictions']:
        assert not (args.output / folder).exists()
    obj = {'hashes': hashes(), 'candidates': candidates(), 'rules': 'DETERMINISTIC_REPAIR.md',
           'encoding_freeze': sha(args.output / 'deterministic_encoding.json'),
           'core_protocol': sha(args.output / 'protocol.json'), 'before_final_test_evaluation': True}
    path = args.output / 'deterministic_protocol.json'
    if path.exists(): assert read(path) == obj
    else: dump(path, obj)
    print('DET_FROZEN', sha(path), flush=True)


def smoke(args, d):
    outer = np.flatnonzero(d.parts != 'external'); tr, va = split_inner(d, outer, 0); rows = []
    for c in candidates()[:4]:
        c = {**c, 'cap': 2, 'schedule': 'cosine'}
        r, cp, pred = ext.train(d, c, tr[:128], va[:32], 11, keep=True)
        y = pred(va[:32]); delta = float(np.abs(ext.replay(d, cp, va[:32][::-1])[::-1] - y).max())
        assert delta < 1e-4
        rows.append({'recipe': c['input_recipe'], 'nin': cp['nin'], 'replay_max': delta})
        del cp, pred
    dump(args.output / 'deterministic_smoke.json', {'status': 'pass', 'recipes': rows})
    print(rows, flush=True)


def main():
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['prepare', 'smoke', 'search', 'final', 'fixed_tree'])
    p.add_argument('--release', type=Path, default=Path('RELEASE_ROOT'))
    p.add_argument('--previous', type=Path, default=Path('BASE_INPUTS'))
    p.add_argument('--output', type=Path, default=HERE / 'results'); p.add_argument('--parts', nargs='+')
    args = p.parse_args(); torch.set_num_threads(1)
    torch.backends.cuda.matmul.allow_tf32 = False; torch.backends.cudnn.allow_tf32 = False
    d = CorrectedDataset(args.release, args.previous); globals()[args.stage](args, d)


if __name__ == '__main__': main()
