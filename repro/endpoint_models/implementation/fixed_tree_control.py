"""Exact native-float32 journal tree control, plus training-only size diagnostic.

The neural Dataset exposes float64 arrays for transforms. This control explicitly
restores the original feature-matrix dtype BEFORE imputation as well as fitting.
Use this file, not the unused float64 fixed_tree stage in deterministic_repair.py.
No neural training/selection source or completed fit is changed.
"""
from pathlib import Path
import argparse
import json
import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from core import dump, outer_partitions, sha, split_inner
from deterministic_repair import CorrectedDataset

HERE=Path(__file__).parent
CONFIG={'n_estimators':360,'min_samples_leaf':2,'max_features':.7}


def read(p):return json.loads(Path(p).read_text())


def tree_arrays(d,response,embedding=True):
    blocks=[d.struct]
    if embedding:blocks.append(d.embeddings['molformer'])
    if response:blocks.append(d.response)
    return np.column_stack(blocks).astype(np.float32)


def model(seed):
    return make_pipeline(SimpleImputer(strategy='median',add_indicator=True),
        ExtraTreesRegressor(**CONFIG,random_state=seed,n_jobs=2))


def protocol(args):
    obj={'script_sha256':sha(__file__),'encoding_freeze_sha256':sha(args.output/'deterministic_encoding.json'),
         'tree':CONFIG,'seeds':[11,29,47],'features_dtype':'float32 before imputation, matching original journal pipeline',
         'imputer':{'strategy':'median','add_indicator':True},'sample_weight':'public 1, ARROW 3',
         'inner_selection':False,'before_final_test_evaluation':True,
         'size_control':'Same original structure/response tree algorithm refit on training-only size split; no held labels fitted.',
         'reason':'Preserve native feature dtype rather than the neural Dataset float64 transform workspace.'}
    path=args.output/'fixed_tree_control_protocol.json'
    if path.exists():assert read(path)==obj
    else:
        assert not (args.output/'comprehensive_policy.json').exists()
        dump(path,obj)
    return obj


def final(args,d):
    freeze=protocol(args);assert (args.output/'comprehensive_policy.json').exists()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        for response in [True,False]:
            name='molformer_fixed_tree_'+('solvai' if response else 'structure');values=[]
            x=tree_arrays(d,response)
            for seed in [11,29,47]:
                path=args.output/'deterministic_tree_final'/part/f'{name}_{seed}.json'
                if path.exists():
                    r=read(path);assert r['control_protocol']==freeze;values.append(np.array(r['predictions']));continue
                fitted=model(seed)
                fitted.fit(x[outer],d.y[outer],extratreesregressor__sample_weight=d.weights[outer])
                y=fitted.predict(x[test]);values.append(y)
                r={'model':name,'seed':seed,'partition':part,'with_response':response,
                   'train_ids':d.frame.id.iloc[outer].tolist(),'validation_ids':[],
                   'test_ids':d.frame.id.iloc[test].tolist(),'predictions':y.tolist(),'control_protocol':freeze}
                if part=='external':
                    cp=args.output/'deterministic_tree_checkpoints'/f'{name}_{seed}.joblib'
                    cp.parent.mkdir(exist_ok=True);joblib.dump(fitted,cp,compress=3)
                    delta=float(np.abs(joblib.load(cp).predict(x[test])-y).max());assert delta<1e-10
                    r.update(checkpoint_sha256=sha(cp),replay_max_error=delta)
                dump(path,r);print('NATIVE_TREE_REFIT',part,name,seed,flush=True)
            frame=d.frame.iloc[test].copy();frame['model']=name;frame['prediction']=np.mean(values,axis=0)
            for seed,y in zip([11,29,47],values):frame[f'seed_{seed}']=y
            path=args.output/'deterministic_tree_predictions'/f'{name}_{part}.csv'
            path.parent.mkdir(exist_ok=True);frame.to_csv(path,index=False)


def size(args,d):
    freeze=protocol(args);assert (args.output/'comprehensive_policy.json').exists()
    outer=np.flatnonzero(d.parts!='external');tr,va=split_inner(d,outer,3)
    for response in [True,False]:
        name='tree_'+('solvai' if response else 'structure')
        path=args.output/'tree_size_diagnostic'/f'{name}.json'
        if path.exists():assert read(path)['control_protocol']==freeze;continue
        x=tree_arrays(d,response,embedding=False);values=[]
        for seed in [11,29,47]:
            fitted=model(seed)
            fitted.fit(x[tr],d.y[tr],extratreesregressor__sample_weight=d.weights[tr])
            values.append(fitted.predict(x[va]));print('SIZE_TREE_REFIT',name,seed,flush=True)
        dump(path,{'model':name,'train_ids':d.frame.id.iloc[tr].tolist(),'validation_ids':d.frame.id.iloc[va].tolist(),
                   'predictions':np.mean(values,axis=0).tolist(),'seed_predictions':[v.tolist() for v in values],
                   'control_protocol':freeze,'selection_used_size_holdout':False})


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','final','size'])
    p.add_argument('--output',type=Path,default=HERE/'results');p.add_argument('--parts',nargs='+')
    args=p.parse_args()
    if args.stage=='prepare':print(protocol(args));return
    d=CorrectedDataset('RELEASE_ROOT',
        'BASE_INPUTS')
    globals()[args.stage](args,d)


if __name__=='__main__':main()
