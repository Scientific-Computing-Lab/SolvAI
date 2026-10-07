"""Replay, homologue stress and exploratory size checks for the capacity option."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import torch
from core import Dataset,SEEDS,fit,dump,sha,split_inner,checkpoint_predict
import tabm_extension as ext
from refine_representation import RefinedDataset
from deterministic_repair import CorrectedDataset
from evaluate import load_stress

HERE=Path(__file__).parent;OUT=HERE/'results'
def read(p):return json.loads(Path(p).read_text())

def get_dataset(args,lock):
    which=lock['choice'].get('dataset','base')
    cls={'base':Dataset,'refined':RefinedDataset,'corrected':CorrectedDataset}[which]
    return cls(args.release,args.previous)

def replay(args,d,lock):
    checks=[]
    if not lock['fallback']:
        ids=np.flatnonzero(d.parts=='external')
        fn=checkpoint_predict if lock['choice']['pipeline']=='core' else ext.replay
        for condition in ['solvai','structure']:
            name='capacity_adapted_'+condition
            for seed in SEEDS:
                path=OUT/'capacity_checkpoints'/f'{name}_{seed}.pt'
                cp=torch.load(path,map_location='cpu',weights_only=False)
                r=read(OUT/'capacity_final/external'/f'{name}_{seed}.json')
                assert cp['config']==lock['choice']['candidate']
                p=fn(d,cp,ids);delta=float(np.abs(p-np.array(r['predictions'])).max());assert delta<1e-4
                checks.append({'model':name,'seed':seed,'delta':delta,'sha256':sha(path)})
        assert len(checks)==10
    dump(OUT/'capacity_checkpoint_replay.json',{'status':'pass','fallback':lock['fallback'],'checks':checks})

def stress(args,d,lock):
    checks=[];rows=[]
    if lock['fallback']:
        f=pd.read_csv(OUT/'stress_predictions.csv',float_precision='round_trip')
        for condition in ['solvai','structure']:
            x=f[f.model==lock['fallback_core_family']+'_'+condition].copy()
            assert len(x)==56;x['model']='capacity_adapted_'+condition;rows.append(x)
    else:
        sd,grid=load_stress(args,d);ids=np.arange(56)
        fn=checkpoint_predict if lock['choice']['pipeline']=='core' else ext.replay
        for condition in ['solvai','structure']:
            name='capacity_adapted_'+condition;values=[]
            for seed in SEEDS:
                path=OUT/'capacity_checkpoints'/f'{name}_{seed}.pt'
                cp=torch.load(path,map_location='cpu',weights_only=False)
                assert cp['config']==lock['choice']['candidate']
                p=fn(sd,cp,ids);q=fn(sd,cp,ids[::-1])[::-1]
                delta=float(np.abs(p-q).max());assert delta<1e-4 and np.isfinite(p).all()
                same=p[(grid.n==0)&grid.series.isin(['glycine','alanine'])]
                assert abs(same[0]-same[1])<1e-4
                values.append(p);checks.append({'model':name,'seed':seed,'permutation_delta':delta,'sha256':sha(path)})
            x=grid.copy();x['model']=name;x['prediction']=np.mean(values,axis=0)
            for seed,p in zip(SEEDS,values):x[f'seed_{seed}']=p
            rows.append(x)
        assert len(checks)==10
    pd.concat(rows,ignore_index=True).to_csv(OUT/'capacity_stress_predictions.csv',index=False)
    dump(OUT/'capacity_stress_replay.json',{'status':'pass','fallback':lock['fallback'],'checks':checks})

def size(args,d,lock):
    tr,va=split_inner(d,np.flatnonzero(d.parts!='external'),3)
    for response in [True,False]:
        suffix='solvai' if response else 'structure';name='capacity_adapted_'+suffix
        path=OUT/'capacity_size_diagnostic'/f'{name}.json'
        if path.exists():continue
        if lock['fallback']:
            r=read(OUT/'size_diagnostic'/f'{lock["fallback_core_family"]}_{suffix}.json')
            r.update(model=name,capacity_core_fallback=True)
        else:
            c=lock['choice']['candidate'];schedule=lock['choice']['refit_schedule'];values=[];records=[]
            for seed in [11,29,47]:
                if lock['choice']['pipeline']=='core':
                    r,cp,pred=fit(d,c,tr,np.array([],dtype=int),seed,response,schedule_override=schedule)
                else:r,cp,pred=ext.train(d,c,tr,np.array([],dtype=int),seed,response,schedule)
                values.append(pred(va));records.append(r);del cp,pred
            r={'model':name,'validation_ids':d.frame.id.iloc[va].tolist(),'predictions':np.mean(values,axis=0).tolist(),
               'seed_predictions':[p.tolist() for p in values],'records':records,'capacity_core_fallback':False,
               'selection_used_size_holdout':False,'caveat':'Exploratory size holdout; labels available during random-fold development.'}
        dump(path,r)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['replay','stress','size'])
    ap.add_argument('--release',type=Path,default=Path('RELEASE_ROOT'))
    ap.add_argument('--previous',type=Path,default=HERE.parent/'neural_endpoint_20261006')
    ap.add_argument('--output',type=Path,default=OUT);args=ap.parse_args()
    assert (OUT/'comprehensive_policy.json').exists()
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    lock=read(OUT/'capacity_locked/external.json');d=get_dataset(args,lock)
    globals()[args.stage](args,d,lock)
