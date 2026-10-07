"""Inner-only capacity boundary challenge; all original searches stay frozen."""
from pathlib import Path
import argparse
import json
import time
import traceback
import numpy as np
import torch
from core import Dataset,SEEDS,fit,dump,sha,outer_partitions,split_inner,schedule_from_records,checkpoint_predict
import tabm_extension as ext
from refine_representation import RefinedDataset
from deterministic_repair import CorrectedDataset,hashes as deterministic_hashes

HERE=Path(__file__).parent
OUT=HERE/'results'
PARTS=[f'arrow_{i}' for i in range(5)]+['external']
FAMILIES=['resnet','dual','ridge_residual','embedding','graph']
EXTRAS={'tabm':('tabm_locked','extension_inner','base'),
        'embedding_refined':('representation_locked','representation_inner','refined'),
        'molformer_deterministic':('deterministic_locked','deterministic_inner','corrected')}

def read(p):return json.loads(Path(p).read_text())
def eligible(c):return c.get('encoder')!='molformer' or c.get('encoder_pipeline')=='deterministic'
def hashes():return {**deterministic_hashes(),**{n:sha(HERE/n) for n in ['capacity_extension.py','CAPACITY_BOUNDARY.md','core.py']}}

def prepare(args,d):
    assert not (OUT/'comprehensive_policy.json').exists()
    arch={}
    for part in PARTS:
        records=read(OUT/'selection'/f'{part}_refine.json')
        arch[part]={f:min([r for r in records if r['candidate']['family']==f],
                       key=lambda r:(r['mean_score'],r['candidate']['id']))['candidate'] for f in FAMILIES}
        for f,(folder,_,_) in EXTRAS.items():arch[part][f]=read(OUT/folder/f'{part}.json')['candidate']
    flagged=[f for f in FAMILIES+list(EXTRAS) if sum(arch[p][f]['width']>=256 for p in PARTS)>=4]
    record={'hashes':hashes(),'core_protocol_sha256':sha(OUT/'protocol.json'),
        'architectures':arch,'flagged_families':flagged,'rules':'CAPACITY_BOUNDARY.md',
        'before_final_evaluation':True,'final_seeds':SEEDS}
    path=OUT/'capacity_protocol.json'
    if path.exists():assert read(path)==record
    else:dump(path,record)
    print('CAPACITY_FROZEN',flagged,sha(path),flush=True)

def inner(d,part,outer,c,j,s,pipeline):
    path=OUT/'capacity_inner'/part/f'{c["id"]}_split{j}_seed{s}.json'
    if path.exists():
        r=read(path);assert r['config']==c and not r.get('failed');return r
    tr,va=split_inner(d,outer,j)
    print('CAPACITY_START',part,c['family'],j,s,flush=True)
    try:
        r,cp,pred=fit(d,c,tr,va,s) if pipeline=='core' else ext.train(d,c,tr,va,s)
        r.update(partition=part,split=j,pipeline=pipeline,protocol_sha256=sha(OUT/'capacity_protocol.json'))
        dump(path,r);del cp,pred
    except Exception:
        dump(path,{'config':c,'failed':True,'traceback':traceback.format_exc()});raise
    print('CAPACITY_DONE',part,c['family'],j,s,r['best_score'],flush=True)
    return r

def search(args,d):
    protocol=read(OUT/'capacity_protocol.json')
    datasets={'base':d,'refined':RefinedDataset(d.release,d.previous),'corrected':CorrectedDataset(d.release,d.previous)}
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        path=OUT/'capacity_locked'/f'{part}.json'
        if path.exists():continue
        lockpath=OUT/'locked'/f'{part}.json'
        while not lockpath.exists():time.sleep(30)
        core=read(lockpath);accepted=[];trials=[]
        groups=[(w['family'],'inner','base','core',w) for w in core['winners']]
        groups += [(f,fitfolder,data,'extension',read(OUT/lockfolder/f'{part}.json'))
                   for f,(lockfolder,fitfolder,data) in EXTRAS.items()]
        for family,fitfolder,data,pipeline,w in groups:
            base=w['candidate']
            if family not in protocol['flagged_families'] or base['width']<256:continue
            if not eligible(base):
                trials.append({'family':family,'base_candidate':base,'retained':False,
                    'excluded_reason':'Invalid legacy stochastic MoLFormer encoding; corrected pipeline has its own comparison.'})
                continue
            c={**base,'id':'wide_'+family,'width':512}
            if family=='graph':c['message_dim']=512
            larger=[inner(datasets[data],part,outer,c,j,47,pipeline) for j in range(3)]
            tag='long' if w['use_long'] else 'ordinary'
            small=[read(OUT/fitfolder/part/f'{base["id"]}_{tag}_split{j}_seed47.json') for j in range(3)]
            diff=np.array([a['best_score']-b['best_score'] for a,b in zip(larger,small)])
            keep=bool(diff.mean()<-.005 and np.count_nonzero(diff<0)>=2)
            trials.append({'family':family,'base_candidate':base,'larger_candidate':c,
                'large_minus_base':diff.tolist(),'retained':keep,'pipeline':pipeline,'dataset':data,'fit_folder':fitfolder})
            if keep:
                rows=[inner(datasets[data],part,outer,c,j,s,pipeline) for j in range(3) for s in [11,29,47]]
                accepted.append({'family':family,'candidate':c,'mean_score':float(np.mean([r['best_score'] for r in rows])),
                    'refit_schedule':schedule_from_records(rows),'chosen_epochs':[r['best_epoch'] for r in rows],
                    'pipeline':pipeline,'dataset':data})
        fallback=min([w for w in core['winners'] if eligible(w['candidate'])],key=lambda w:(w['mean_score'],w['family']))
        choice=min(accepted,key=lambda w:(w['mean_score'],w['family'])) if accepted else fallback
        dump(path,{'partition':part,'choice':choice,'fallback':not bool(accepted),
            'fallback_core_family':fallback['family'],'trials':trials,'retained':accepted,
            'protocol_sha256':sha(OUT/'capacity_protocol.json'),'test_labels_used':False})
        print('CAPACITY_LOCKED',part,choice['family'],bool(accepted),flush=True)

def final(args,d):
    assert (OUT/'comprehensive_policy.json').exists()
    datasets={'base':d,'refined':RefinedDataset(d.release,d.previous),'corrected':CorrectedDataset(d.release,d.previous)}
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        w=read(OUT/'capacity_locked'/f'{part}.json')
        if w['fallback']:continue
        c=w['choice']['candidate'];schedule=w['choice']['refit_schedule']
        dataset=datasets[w['choice']['dataset']];pipeline=w['choice']['pipeline']
        for response in [True,False]:
            name='capacity_adapted_'+('solvai' if response else 'structure');values=[]
            for seed in SEEDS:
                path=OUT/'capacity_final'/part/f'{name}_{seed}.json'
                if path.exists():values.append(np.array(read(path)['predictions']));continue
                if pipeline=='core':
                    r,cp,pred=fit(dataset,c,outer,np.array([],dtype=int),seed,response,schedule_override=schedule,checkpoint=part=='external')
                else:r,cp,pred=ext.train(dataset,c,outer,np.array([],dtype=int),seed,response,schedule,part=='external')
                p=pred(test);assert np.isfinite(p).all()
                r.update(model=name,partition=part,test_ids=d.frame.id.iloc[test].tolist(),predictions=p.tolist(),
                         pipeline=pipeline,dataset=w['choice']['dataset'],protocol_sha256=sha(OUT/'capacity_protocol.json'))
                if cp is not None:
                    cpp=OUT/'capacity_checkpoints'/f'{name}_{seed}.pt';cpp.parent.mkdir(exist_ok=True)
                    torch.save(cp,cpp);q=checkpoint_predict(dataset,cp,test) if pipeline=='core' else ext.replay(dataset,cp,test)
                    delta=float(np.abs(q-p).max());assert delta<1e-4
                    r.update(checkpoint_sha256=sha(cpp),replay_max_error=delta)
                dump(path,r);values.append(p);del cp,pred
            f=d.frame.iloc[test].copy();f['model']=name;f['prediction']=np.mean(values,axis=0)
            for seed,p in zip(SEEDS,values):f[f'seed_{seed}']=p
            path=OUT/'capacity_predictions'/f'{name}_{part}.csv';path.parent.mkdir(exist_ok=True)
            f.to_csv(path,index=False)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','search','final'])
    ap.add_argument('--parts',nargs='+');args=ap.parse_args()
    torch.set_num_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    d=Dataset('RELEASE_ROOT',HERE.parent/'neural_endpoint_20261006')
    if args.stage!='prepare':assert read(OUT/'capacity_protocol.json')['hashes']==hashes()
    globals()[args.stage](args,d)

if __name__=='__main__':main()
