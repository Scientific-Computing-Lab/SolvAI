"""Pinned, offline TabPFN extension; no changes to any released model or data."""
from pathlib import Path
import os
import sys
HERE=Path(__file__).parent
sys.path.insert(0,str(HERE/'tabpfn_vendor'))
os.environ['HF_HUB_OFFLINE']='1'
os.environ['TRANSFORMERS_OFFLINE']='1'
import argparse
import gc
import hashlib
import importlib.metadata
import json
import time
import traceback
import numpy as np
import pandas as pd
import torch
from tabpfn import TabPFNRegressor
from tabpfn.model_loading import save_fitted_tabpfn_model, load_fitted_tabpfn_model
from core import SEEDS, dump, outer_partitions, sha, split_inner, scores
from deterministic_repair import CorrectedDataset

DATA=HERE/'results'
PARTS=[f'arrow_{i}' for i in range(5)]+['external']
CHECKPOINTS={'full':'tabpfn-v3.5-20260909.safetensors','fast':'tabpfn-v3.5-fast-20260909.safetensors'}

def read(p):return json.loads(Path(p).read_text())

def candidates():
    return [{'id':f'p{i:02d}','checkpoint':ckpt,'recipe':recipe,'n_estimators':8,'arrow_repeats':3}
            for i,(ckpt,recipe) in enumerate((a,b) for a in ['full','fast']
                for b in ['binary','raw_count','no_fingerprint','molformer_no_fingerprint'])]

def hashes():
    paths=[HERE/n for n in ['tabpfn_extension.py','TABPFN_EXTENSION.md','core.py',
           'tabm_extension.py','refine_representation.py','deterministic_repair.py',
           'results/deterministic_encoding.json','results/tabpfn_assets.json']]
    paths+=sorted((HERE/'tabpfn_vendor').rglob('*.py'))
    paths+=sorted((HERE/'tabpfn_vendor').glob('*.dist-info/METADATA'))
    return {str(p.relative_to(HERE)):sha(p) for p in paths}

def no_network(event,args):
    if event in ['socket.connect','socket.getaddrinfo']:
        raise RuntimeError('Network disabled during local TabPFN research inference')

def arrays(d,c,response=True):
    recipe=c['recipe'];embedded=recipe=='molformer_no_fingerprint'
    cfg={'family':'embedding' if embedded else 'resnet','encoder':'molformer',
         'representation':'binary_log' if recipe=='binary' else 'count_log',
         'input_recipe':'no_fingerprint' if 'no_fingerprint' in recipe else recipe}
    x,_,_=d.arrays(cfg,response)
    # Leave missing values to the official, training-fitted preprocessing.
    return np.where(np.isfinite(x),x,np.nan).astype(np.float32)

def fitted_mean_support(model):
    # FullSupportBarDistribution.mean is a convex combination of these finite
    # bucket means, even though the edge-bin distributions have unbounded tails.
    b=model.znorm_space_bardist_;means=b.borders[:-1]+b.bucket_widths/2
    means=means.clone()
    means[0]=b.borders[1]-b.halfnormal_with_p_weight_before(b.bucket_widths[0]).mean
    means[-1]=b.borders[-2]+b.halfnormal_with_p_weight_before(b.bucket_widths[-1]).mean
    means=means*model.y_train_std_+model.y_train_mean_
    return {'lower':float(means.min().cpu()),'upper':float(means.max().cpu()),
            'interpretation':'Finite convex-hull bound for mean output, not quantile support and not necessarily tight or attainable.'}

def fit_model(d,c,tr,seed,response=True):
    assert not np.any(d.parts[tr]=='external')
    x=arrays(d,c,response)
    multiplicities=np.where(d.parts[tr]=='arrow',c['arrow_repeats'],1)
    context=np.repeat(tr,multiplicities)
    assert set(context)==set(tr)
    model=TabPFNRegressor(model_path=str(HERE/'tabpfn_assets'/CHECKPOINTS[c['checkpoint']]),
        n_estimators=c['n_estimators'],device='cuda',fit_mode='fit_preprocessors',
        random_state=seed,n_preprocessing_jobs=1,show_progress_bar=False)
    start=time.monotonic();model.fit(x[context],d.y[context])
    record={'config':c,'seed':seed,'with_response':response,'train_ids':d.frame.id.iloc[tr].tolist(),
        'context_rows':len(context),'unique_training_rows':len(tr),'fit_seconds':time.monotonic()-start,
        'input_features':x.shape[1],'n_estimators_actual':model.n_estimators_,'validation_ids':[],
        'mean_output_support':fitted_mean_support(model),'inference_config':str(model.inference_config_),
        'protocol_sha256':sha(DATA/'tabpfn_protocol.json')}
    return model,x,record

def clear(model):
    del model
    gc.collect();torch.cuda.empty_cache()

def inner(d,part,outer,c,split,seed):
    name=f'{c["id"]}_r{c["arrow_repeats"]}_e{c["n_estimators"]}_split{split}_seed{seed}.json'
    path=DATA/'tabpfn_inner'/part/name
    if path.exists():
        r=read(path);assert r['config']==c and not r.get('failed');return r
    tr,va=split_inner(d,outer,split)
    print('PFN_START',part,name,flush=True)
    try:
        model,x,r=fit_model(d,c,tr,seed)
        start=time.monotonic();p=model.predict(x[va]).astype(float)
        assert np.isfinite(p).all()
        r.update(partition=part,split=split,validation_ids=d.frame.id.iloc[va].tolist(),
                 validation_predictions=p.tolist(),best_score=scores(d,va,p)['score'],
                 validation_scores=scores(d,va,p),predict_seconds=time.monotonic()-start)
        dump(path,r);del model;gc.collect();torch.cuda.empty_cache()
        print('PFN_DONE',part,name,round(r['best_score'],6),flush=True)
        return r
    except Exception:
        dump(path,{'config':c,'failed':True,'traceback':traceback.format_exc()});raise

def improve(diffs):return bool(np.mean(diffs)<-.005 and np.count_nonzero(np.array(diffs)<0)>=2)

def search(args,d):
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        path=DATA/'tabpfn_locked'/f'{part}.json'
        if path.exists():continue
        screen=[(inner(d,part,outer,c,0,11)['best_score'],c) for c in candidates()]
        repeated=[]
        for _,c in sorted(screen,key=lambda a:(a[0],a[1]['id']))[:2]:
            r=[inner(d,part,outer,c,j,s) for j in range(3) for s in [11,29]]
            repeated.append((float(np.mean([a['best_score'] for a in r])),c))
        _,base=min(repeated,key=lambda a:(a[0],a[1]['id']))
        uniform={**base,'arrow_repeats':1}
        normal=[inner(d,part,outer,base,j,47) for j in range(3)]
        unweighted=[inner(d,part,outer,uniform,j,47) for j in range(3)]
        changes=[b['best_score']-a['best_score'] for a,b in zip(normal,unweighted)]
        chosen=uniform if improve(changes) else base
        baseline=[inner(d,part,outer,chosen,j,s) for j in range(3) for s in [11,29,47]]
        boosted={**chosen,'n_estimators':32}
        bigger=[inner(d,part,outer,boosted,j,47) for j in range(3)]
        changes_compute=[b['best_score']-inner(d,part,outer,chosen,j,47)['best_score'] for j,b in enumerate(bigger)]
        if improve(changes_compute):
            chosen=boosted
            baseline=[inner(d,part,outer,chosen,j,s) for j in range(3) for s in [11,29,47]]
        dump(path,{'candidate':chosen,'base_candidate':base,
            'mean_score':float(np.mean([r['best_score'] for r in baseline])),
            'uniform_minus_repeated':changes,'boosted_minus_eight':changes_compute,
            'uniform_selected':improve(changes),'boosted_selected':improve(changes_compute),
            'protocol_sha256':sha(DATA/'tabpfn_protocol.json'),'test_labels_used':False})
        print('PFN_LOCKED',part,chosen,flush=True)

def final(args,d):
    assert (DATA/'comprehensive_policy.json').exists()
    for part in PARTS:assert (DATA/'tabpfn_locked'/f'{part}.json').exists()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        c=read(DATA/'tabpfn_locked'/f'{part}.json')['candidate']
        for response in [True,False]:
            name='tabpfn_'+('solvai' if response else 'structure');values=[]
            for seed in SEEDS:
                path=DATA/'tabpfn_final'/part/f'{name}_{seed}.json'
                if path.exists():values.append(np.array(read(path)['predictions']));continue
                model,x,r=fit_model(d,c,outer,seed,response)
                start=time.monotonic();p=model.predict(x[test]).astype(float)
                assert np.isfinite(p).all()
                r.update(model=name,partition=part,test_ids=d.frame.id.iloc[test].tolist(),
                    validation_ids=[],predictions=p.tolist(),predict_seconds=time.monotonic()-start)
                if part=='external':
                    cp=DATA/'tabpfn_private_context'/f'{name}_{seed}.tabpfn_fit'
                    save_fitted_tabpfn_model(model,cp)
                    r['checkpoint_sha256']=sha(cp)
                    del model;gc.collect();torch.cuda.empty_cache()
                    model=load_fitted_tabpfn_model(cp,device='cuda')
                    q=model.predict(x[test]).astype(float)
                    delta=float(np.abs(p-q).max());assert delta<1e-4,delta
                    r['replay_max_error']=delta
                dump(path,r);values.append(p);del model;gc.collect();torch.cuda.empty_cache()
                print('PFN_FINAL',part,name,seed,flush=True)
            f=d.frame.iloc[test].copy();f['model']=name;f['prediction']=np.mean(values,axis=0)
            for seed,p in zip(SEEDS,values):f[f'seed_{seed}']=p
            path=DATA/'tabpfn_predictions'/f'{name}_{part}.csv';path.parent.mkdir(exist_ok=True)
            f.to_csv(path,index=False)

def stress(args,d):
    from evaluate import load_stress
    sd,grid=load_stress(args,d)
    c=read(DATA/'tabpfn_locked/external.json')['candidate'];rows=[];checks=[]
    for response in [True,False]:
        name='tabpfn_'+('solvai' if response else 'structure');x=arrays(sd,c,response);values=[]
        for seed in SEEDS:
            cp=DATA/'tabpfn_private_context'/f'{name}_{seed}.tabpfn_fit'
            model=load_fitted_tabpfn_model(cp,device='cuda')
            p=model.predict(x).astype(float);q=model.predict(x[::-1]).astype(float)[::-1]
            delta=float(np.abs(p-q).max());assert delta<1e-3,delta
            same=p[(grid.n==0)&grid.series.isin(['glycine','alanine'])]
            assert abs(same[0]-same[1])<1e-4
            assert np.isfinite(p).all()
            values.append(p);checks.append({'model':name,'seed':seed,'permutation_max_error':delta,
                'mean_output_support':fitted_mean_support(model),'checkpoint_sha256':sha(cp)})
            del model;gc.collect();torch.cuda.empty_cache()
        f=grid.copy();f['model']=name;f['prediction']=np.mean(values,axis=0)
        for seed,p in zip(SEEDS,values):f[f'seed_{seed}']=p
        rows.append(f)
    pd.concat(rows,ignore_index=True).to_csv(DATA/'tabpfn_stress_predictions.csv',index=False)
    dump(DATA/'tabpfn_replay.json',{'status':'pass','checks':checks,'scope':'No peptide truth or calibrated uncertainty.'})

def size(args,d):
    c=read(DATA/'tabpfn_locked/external.json')['candidate']
    tr,va=split_inner(d,np.flatnonzero(d.parts!='external'),3)
    for response in [True,False]:
        name='tabpfn_'+('solvai' if response else 'structure')
        path=DATA/'tabpfn_size_diagnostic'/f'{name}.json'
        if path.exists():continue
        values=[];records=[]
        for seed in [11,29,47]:
            model,x,r=fit_model(d,c,tr,seed,response)
            values.append(model.predict(x[va]).astype(float));records.append(r)
            del model;gc.collect();torch.cuda.empty_cache()
        dump(path,{'model':name,'validation_ids':d.frame.id.iloc[va].tolist(),
            'predictions':np.mean(values,axis=0).tolist(),'seed_predictions':[p.tolist() for p in values],
            'records':records,'selection_used_size_holdout':False,
            'caveat':'Exploratory size holdout; labels available during random-split development.'})

def prepare(args,d):
    assert not (DATA/'comprehensive_policy.json').exists()
    for folder in ['predictions','extension_predictions','representation_predictions','deterministic_predictions',
                   'deterministic_tree_predictions','tabpfn_predictions']:
        assert not (DATA/folder).exists(),folder
    record={'hashes':hashes(),'candidates':candidates(),'tabpfn_version':importlib.metadata.version('tabpfn'),
        'final_seeds':SEEDS,'assets':read(DATA/'tabpfn_assets.json'),'source_weights':d.source_hashes,
        'source_protocol':sha(DATA/'protocol.json'),'selection_metric':'equal-source validation MAE',
        'before_final_test_evaluation':True,'rules':'TABPFN_EXTENSION.md','private_context_stays_on_clark':True}
    path=DATA/'tabpfn_protocol.json'
    if path.exists():assert read(path)==record
    else:dump(path,record)
    print('PFN_FROZEN',sha(path),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','search','final','stress','size'])
    ap.add_argument('--parts',nargs='+')
    ap.add_argument('--release',type=Path,default=Path('RELEASE_ROOT'))
    ap.add_argument('--previous',type=Path,default=HERE.parent/'neural_endpoint_20261006')
    ap.add_argument('--output',type=Path,default=DATA)
    args=ap.parse_args()
    torch.set_num_threads(1);torch.cuda.set_per_process_memory_fraction(.4)
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    for row in read(DATA/'tabpfn_assets.json')['files']:
        assert sha(HERE/'tabpfn_assets'/row['filename'])==row['sha256']
    d=CorrectedDataset(args.release,args.previous)
    sys.addaudithook(no_network)
    if args.stage!='prepare':assert read(DATA/'tabpfn_protocol.json')['hashes']==hashes()
    if args.stage in ['stress','size']:assert (DATA/'comprehensive_policy.json').exists()
    globals()[args.stage](args,d)

if __name__=='__main__':main()
