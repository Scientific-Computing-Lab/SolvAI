"""Inner-validation-driven input refinement; retains both preceding studies."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import torch

import tabm_extension as ext
from core import Dataset, SEEDS, dump, outer_partitions, schedule_from_records, sha, split_inner

HERE=Path(__file__).parent


def read(p):return json.loads(Path(p).read_text())


class RefinedDataset(Dataset):
    def arrays(self,c,response=True):
        x,b,binary=super().arrays(c,response)
        if c.get('input_recipe')=='no_fingerprint':
            keep=np.r_[np.arange(217),np.arange(2265,x.shape[1])]
            x=x[:,keep];binary=binary[keep];b-=2048
        return x,b,binary


def hashes():
    return {**ext.hashes(),'refine_representation.py':sha(HERE/'refine_representation.py'),
            'REPRESENTATION_REFINEMENT.md':sha(HERE/'REPRESENTATION_REFINEMENT.md')}


def candidates(args,part):
    original=read(args.output/'selection'/f'{part}_screen.json');out=[]
    for encoder in ['molformer','unimol']:
        source=min([r for r in original if r['candidate'].get('encoder')==encoder],
                   key=lambda r:(r['score'],r['candidate']['id']))['candidate']
        for recipe in ['binary','raw_count','log_count','no_fingerprint']:
            out.append({**source,'id':f'e{len(out):02d}','representation':'binary_log' if recipe=='binary' else 'count_log',
                        'count_scaling':recipe,'input_recipe':recipe,'source_config_id':source['id']})
    return out


def inner(args,d,part,outer,c,split,seed,tag='ordinary'):
    p=args.output/'representation_inner'/part/f'{c["id"]}_{tag}_split{split}_seed{seed}.json'
    if p.exists():
        r=read(p);assert r['config']==c and r['refinement_hashes']==hashes();return r
    print('REP_START',part,c['id'],c['encoder'],c['input_recipe'],tag,split,seed,flush=True)
    tr,va=split_inner(d,outer,split);r,cp,pred=ext.train(d,c,tr,va,seed)
    r.update(partition=part,split=split,tag=tag,refinement_hashes=hashes());dump(p,r)
    del cp,pred
    print('REP_DONE',part,c['id'],tag,split,seed,round(r['best_score'],5),r['best_epoch'],flush=True)
    return r


def search(args,d):
    assert read(args.output/'representation_protocol.json')['hashes']==hashes()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        if (args.output/'representation_locked'/f'{part}.json').exists():continue
        path=args.output/'selection'/f'{part}_screen.json'
        while not path.exists():
            print('REP_WAITING_FOR_CORE_SCREEN',part,flush=True);time.sleep(30)
        configs=candidates(args,part);dump(args.output/'representation_candidates'/f'{part}.json',configs)
        screen=[]
        for c in configs:screen.append((inner(args,d,part,outer,c,0,11)['best_score'],c))
        refined=[]
        for _,c in sorted(screen,key=lambda r:(r[0],r[1]['id']))[:2]:
            rr=[inner(args,d,part,outer,c,j,s) for j in range(3) for s in [11,29]]
            refined.append((float(np.mean([r['best_score'] for r in rr])),c))
        _,c=min(refined,key=lambda a:(a[0],a[1]['id']));longer={**c,'cap':1200,'schedule':'cosine'}
        normal=[inner(args,d,part,outer,c,j,47) for j in range(3)]
        long=[inner(args,d,part,outer,longer,j,47,'long') for j in range(3)]
        diffs=np.array([b['best_score']-a['best_score'] for a,b in zip(normal,long)])
        use_long=bool(diffs.mean()<-.005 and (diffs<0).sum()>=2)
        if use_long:
            c=longer;records=long+[inner(args,d,part,outer,c,j,s,'long') for j in range(3) for s in [11,29]]
        else:records=[inner(args,d,part,outer,c,j,s) for j in range(3) for s in [11,29,47]]
        dump(args.output/'representation_locked'/f'{part}.json',{'candidate':c,'mean_score':float(np.mean([r['best_score'] for r in records])),
            'use_long':use_long,'long_minus_normal':diffs.tolist(),'refit_schedule':schedule_from_records(records),
            'chosen_epochs':[r['best_epoch'] for r in records],'refinement_hashes':hashes(),'test_labels_used':False})
        print('REP_LOCKED',part,c['id'],c['input_recipe'],flush=True)


def final(args,d):
    for part,_,_ in outer_partitions(d):
        for folder in ['locked','extension_locked','representation_locked']:assert (args.output/folder/f'{part}.json').exists()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        choice=read(args.output/'representation_locked'/f'{part}.json')
        for response in [True,False]:
            name='embedding_refined_'+('solvai' if response else 'structure');values=[]
            for seed in SEEDS:
                path=args.output/'representation_final'/part/f'{name}_{seed}.json'
                if path.exists():values.append(np.array(read(path)['predictions']));continue
                r,cp,pred=ext.train(d,choice['candidate'],outer,np.array([],dtype=int),seed,response,choice['refit_schedule'],part=='external')
                y=pred(test);r.update(model=name,partition=part,test_ids=d.frame.id.iloc[test].tolist(),predictions=y.tolist(),refinement_hashes=hashes())
                if cp is not None:
                    cp_path=args.output/'representation_checkpoints'/f'{name}_{seed}.pt';cp_path.parent.mkdir(exist_ok=True);torch.save(cp,cp_path)
                    delta=float(np.max(np.abs(ext.replay(d,cp,test)-y)));assert delta<1e-4
                    r.update(checkpoint_sha256=sha(cp_path),replay_max_error=delta)
                dump(path,r);values.append(y);del cp,pred
                print('REP_REFIT',part,name,seed,flush=True)
            x=d.frame.iloc[test].copy();x['model']=name;x['prediction']=np.mean(values,axis=0)
            for seed,y in zip(SEEDS,values):x[f'seed_{seed}']=y
            p=args.output/'representation_predictions'/f'{name}_{part}.csv';p.parent.mkdir(exist_ok=True);x.to_csv(p,index=False)


def prepare(args,d):
    assert not (args.output/'predictions').exists()
    p=args.output/'representation_protocol.json'
    obj={'hashes':hashes(),'rules':'see REPRESENTATION_REFINEMENT.md','core_protocol':sha(args.output/'protocol.json'),
         'extension_protocol':sha(args.output/'extension_protocol.json'),'adaptive_from_inner_validation':True,
         'before_final_test_evaluation':True,'recipes':['binary','raw_count','log_count','no_fingerprint']}
    if p.exists():assert read(p)==obj
    else:dump(p,obj)
    print('REP_FROZEN',sha(p),flush=True)


def smoke(args,d):
    from core import candidates as originals
    outer=np.flatnonzero(d.parts!='external');tr,va=split_inner(d,outer,0);out=[]
    for recipe in ['binary','raw_count','log_count','no_fingerprint']:
        c={**originals()[36],'id':'smoke_'+recipe,'input_recipe':recipe,'count_scaling':recipe,
           'representation':'binary_log' if recipe=='binary' else 'count_log','cap':2,'schedule':'cosine'}
        x,b,m=d.arrays(c);assert x.shape[1]==(1000 if recipe=='no_fingerprint' else 3048)
        r,cp,pred=ext.train(d,c,tr[:128],va[:32],11,keep=True)
        a=pred(va[:32]);z=ext.replay(d,cp,va[:32][::-1])[::-1];assert np.max(np.abs(a-z))<1e-4
        out.append({'recipe':recipe,'columns':x.shape[1],'replay_and_permutation':True});del cp,pred
    dump(args.output/'representation_smoke.json',out);print(out,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','smoke','search','final'])
    p.add_argument('--release',type=Path,default=Path('RELEASE_ROOT'))
    p.add_argument('--previous',type=Path,default=Path('BASE_INPUTS'))
    p.add_argument('--output',type=Path,default=HERE/'results');p.add_argument('--parts',nargs='+');args=p.parse_args()
    torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    d=RefinedDataset(args.release,args.previous);globals()[args.stage](args,d)


if __name__=='__main__':main()
