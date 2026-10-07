from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import torch
from torch import nn

HERE=Path(__file__).parent
sys.path.insert(0,str(HERE/'vendor'))
from tabm import TabM
from rtdl_num_embeddings import LinearReLUEmbeddings
from core import (Dataset, Head, Preprocessor, SEEDS, dump, outer_partitions,
                  schedule_from_records, scores, seed_all, sha, split_inner)


def read(p):return json.loads(Path(p).read_text())


def configs():
    rows=[]
    settings=[
        dict(width=128,depth=2,dropout=0.,lr=.002,decay=.001,embedding=0,schedule='plateau'),
        dict(width=256,depth=3,dropout=.1,lr=.001,decay=.003,embedding=0,schedule='plateau'),
        dict(width=64,depth=2,dropout=.1,lr=.001,decay=.001,embedding=4,schedule='cosine'),
        dict(width=128,depth=3,dropout=.1,lr=.0005,decay=.003,embedding=4,schedule='cosine')]
    for mode in ['binary','raw_count','log_count']:
        for setting in settings:
            rows.append(dict(family='tabm',id=f't{len(rows):02d}',representation='binary_log' if mode=='binary' else 'count_log',
                             count_scaling=mode,cap=600,k=32,batch=128,loss='mse',**setting))
    return rows


def hashes():
    paths=[HERE/'core.py',HERE/'tabm_extension.py',HERE/'EXTENSION.md']
    paths += sorted((HERE/'vendor').glob('*.py'))
    paths += sorted((HERE/'vendor').glob('*.dist-info/METADATA'))
    return {str(p.relative_to(HERE)):sha(p) for p in paths}


def arrays(d,c,response):
    x,b,binary=d.arrays(c,response)
    if c.get('count_scaling') in ['raw_count','log_count']:
        assert c['family']!='graph'
        x[:,217:2265]=d.counts if c['count_scaling']=='raw_count' else np.log1p(d.counts)
        # The preprocessor's mask denotes columns kept in their native scale;
        # it is also applicable to counts, not a claim that counts are binary.
        binary[217:2265]=True
    return x,b,binary


class Model(nn.Module):
    def __init__(self,c,nin,boundary):
        super().__init__();self.c=c
        if c['family']=='tabm':
            emb=LinearReLUEmbeddings(nin,c['embedding']) if c['embedding'] else None
            self.net=TabM.make(n_num_features=nin,d_out=1,num_embeddings=emb,k=c['k'],
                               d_block=c['width'],n_blocks=c['depth'],dropout=c['dropout'])
        else:self.net=Head(c,nin,boundary)
    def forward(self,x):
        if self.c['family']=='tabm':return self.net(x).squeeze(-1)
        return self.net(x).unsqueeze(1)


def train(d,c,tr,va,seed,response=True,schedule=None,keep=False):
    from sklearn.linear_model import Ridge
    assert not np.any(d.parts[tr]=='external') and not np.any(d.parts[va]=='external')
    assert not np.intersect1d(tr,va).size
    seed_all(seed);start=time.monotonic()
    x,b,binary=arrays(d,c,response);prep=Preprocessor().fit(x[tr],b,binary);xx=prep.transform(x)
    mean=float(d.y[tr].mean());scale=float(d.y[tr].std());yy=(d.y-mean)/scale
    base=np.zeros(len(yy));ridge_state=None
    if c['family']=='ridge_residual':
        cols=prep.ridge_indices
        r=Ridge(alpha=10.,solver='cholesky').fit(xx[tr][:,cols],yy[tr],sample_weight=d.weights[tr])
        base=r.predict(xx[:,cols]).astype(float)
        ridge_state={'columns':cols,'coef':r.coef_,'intercept':float(r.intercept_)}
    model=Model(c,xx.shape[1],prep.boundary).cuda()
    xt=torch.as_tensor(xx,device='cuda');yt=torch.as_tensor((yy-base).astype(np.float32),device='cuda')
    wt=torch.as_tensor((d.weights/d.weights[tr].mean()).astype(np.float32),device='cuda')
    opt=torch.optim.AdamW([{'params':[p for p in model.parameters() if p.ndim>=2],'weight_decay':c['decay']},
                          {'params':[p for p in model.parameters() if p.ndim<2],'weight_decay':0.}],lr=c['lr'])
    plateau=torch.optim.lr_scheduler.ReduceLROnPlateau(opt,factor=.5,patience=25,threshold=0.,min_lr=1e-6)
    def predict(ids):
        model.eval();out=[]
        with torch.inference_mode():
            for batch in np.array_split(ids,max(1,math.ceil(len(ids)/128))):
                out.extend(model(xt[batch]).mean(1).cpu().numpy().astype(float))
        return (np.array(out)+base[ids])*scale+mean
    trace=[];best=float('inf');best_epoch=0;best_pred=None
    cap=len(schedule) if schedule is not None else c['cap'];rng=np.random.default_rng(seed)
    for epoch in range(1,cap+1):
        if schedule is not None:lr=schedule[epoch-1]
        elif c['schedule']=='cosine':lr=c['lr']*min(1.,epoch/10.)*(.01+.99*(1+math.cos(math.pi*(epoch-1)/cap))/2)
        else:lr=opt.param_groups[0]['lr']
        for group in opt.param_groups:group['lr']=lr
        model.train()
        for batch in np.array_split(rng.permutation(tr),math.ceil(len(tr)/c['batch'])):
            opt.zero_grad(set_to_none=True)
            delta=model(xt[batch])-yt[batch,None]
            losses=delta.square() if c['loss']=='mse' else nn.functional.huber_loss(delta,torch.zeros_like(delta),delta=1.,reduction='none')
            # Crucial: average member LOSSES, never square the mean prediction error.
            loss=(losses.mean(1)*wt[batch]).mean()
            assert torch.isfinite(loss)
            loss.backward();nn.utils.clip_grad_norm_(model.parameters(),5.);opt.step()
        if len(va):
            pred=predict(va);sc=scores(d,va,pred);row={'epoch':epoch,'lr':float(lr),**sc}
            if epoch==1 or epoch%10==0:row['train']=scores(d,tr,predict(tr))
            trace.append(row)
            if sc['score']<best:best=sc['score'];best_epoch=epoch;best_pred=pred.copy()
            if epoch%100==0:print('EXT_EPOCH',c['id'],seed,epoch,round(sc['score'],5),flush=True)
            if c['schedule']=='plateau':plateau.step(sc['score'])
            if c['schedule']=='plateau' and epoch>=150 and epoch-best_epoch>=100:break
    rec={'config':c,'seed':seed,'with_response':response,'epochs_run':epoch,'best_epoch':best_epoch if len(va) else epoch,
         'best_score':best if len(va) else None,'trace':trace,'validation_predictions':best_pred.tolist() if best_pred is not None else [],
         'train_ids':d.frame.id.iloc[tr].tolist(),'validation_ids':d.frame.id.iloc[va].tolist(),
         'train_final':scores(d,tr,predict(tr)),'seconds':time.monotonic()-start,'source_hashes':hashes()}
    cp={'config':c,'state':{k:v.detach().cpu() for k,v in model.state_dict().items()},'prep':prep.state(),
        'mean':mean,'scale':scale,'ridge':ridge_state,'nin':xx.shape[1],'boundary':prep.boundary,'response':response,'seed':seed} if keep else None
    return rec,cp,predict


def replay(d,cp,ids):
    x,_,_=arrays(d,cp['config'],cp['response']);p=Preprocessor();p.__dict__.update(cp['prep']);x=p.transform(x)
    m=Model(cp['config'],cp['nin'],cp['boundary']).cuda();m.load_state_dict(cp['state']);m.eval();result=[]
    with torch.inference_mode():
        for batch in np.array_split(ids,max(1,math.ceil(len(ids)/128))):
            z=m(torch.as_tensor(x[batch],device='cuda')).mean(1).cpu().numpy().astype(float)
            if cp['ridge'] is not None:
                r=cp['ridge'];z+=x[batch][:,r['columns']]@r['coef']+r['intercept']
            result.extend(z*cp['scale']+cp['mean'])
    return np.array(result)


def inner(args,d,part,outer,c,split,seed,tag='ordinary'):
    path=args.output/'extension_inner'/part/f'{c["id"]}_{tag}_split{split}_seed{seed}.json'
    if path.exists():
        r=read(path);assert r['source_hashes']==hashes() and r['config']==c;return r
    print('EXT_START',part,c['id'],tag,split,seed,flush=True)
    tr,va=split_inner(d,outer,split);r,cp,pred=train(d,c,tr,va,seed,schedule=c.get('fixed_schedule'))
    del cp,pred
    r.update(partition=part,split=split,tag=tag);dump(path,r)
    print('EXT_DONE',part,c['id'],tag,split,seed,round(r['best_score'],5),r['best_epoch'],round(r['seconds'],1),flush=True)
    return r


def search(args,d):
    assert read(args.output/'extension_protocol.json')['hashes']==hashes()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        locked=args.output/'extension_locked'/f'{part}.json'
        if locked.exists():continue
        screen=[]
        for c in configs():screen.append((inner(args,d,part,outer,c,0,11)['best_score'],c))
        refined=[]
        for _,c in sorted(screen,key=lambda a:(a[0],a[1]['id']))[:2]:
            rows=[inner(args,d,part,outer,c,j,s) for j in range(3) for s in [11,29]]
            refined.append((np.mean([r['best_score'] for r in rows]),c))
        _,c=min(refined,key=lambda a:(a[0],a[1]['id']));longer={**c,'cap':1200,'schedule':'cosine'}
        normal=[inner(args,d,part,outer,c,j,47) for j in range(3)]
        long=[inner(args,d,part,outer,longer,j,47,'long') for j in range(3)]
        diffs=np.array([b['best_score']-a['best_score'] for a,b in zip(normal,long)])
        use_long=bool(diffs.mean()<-.005 and (diffs<0).sum()>=2)
        if use_long:
            c=longer;rows=long+[inner(args,d,part,outer,c,j,s,'long') for j in range(3) for s in [11,29]]
        else:rows=[inner(args,d,part,outer,c,j,s) for j in range(3) for s in [11,29,47]]
        tabm={'candidate':c,'mean_score':float(np.mean([r['best_score'] for r in rows])),
              'refit_schedule':schedule_from_records(rows),'long_minus_normal':diffs.tolist(),'use_long':use_long}
        dump(args.output/'tabm_locked'/f'{part}.json',tabm)
        # This optional continuation needs the original search's inner-only choice.
        corepath=args.output/'locked'/f'{part}.json'
        if not corepath.exists():
            print('TABM_LOCKED_AWAITING_CORE',part,flush=True);continue
        complete_count_adaptation(args,d,part,outer,tabm)


def complete_count_adaptation(args,d,part,outer,tabm):
    core=read(args.output/'locked'/f'{part}.json')
    base=next(w for w in core['winners'] if w['family']==core['overall_family'])
    adapted=[]
    if base['family']!='graph':
        for mode in ['raw_count','log_count']:
            c={**base['candidate'],'id':'adapt_'+mode,'representation':'count_log','count_scaling':mode,
               'fixed_schedule':base['refit_schedule'],'cap':len(base['refit_schedule'])}
            rows=[inner(args,d,part,outer,c,j,s,'adapt') for j in range(3) for s in [11,29,47]]
            adapted.append({'candidate':c,'mean_score':float(np.mean([r['best_score'] for r in rows])),
                            'refit_schedule':schedule_from_records(rows)})
    adaptation=min(adapted,key=lambda a:a['mean_score']) if adapted else None
    options={'core':base['mean_score'],'tabm':tabm['mean_score']}
    if adaptation is not None:options['count_adapted']=adaptation['mean_score']
    dump(args.output/'extension_locked'/f'{part}.json',{'tabm':tabm,'count_adapted':adaptation,'core_family':base['family'],
         'inner_scores':options,'selected':min(options,key=options.get),'test_scores_used':False,'source_hashes':hashes()})
    print('EXT_LOCKED',part,min(options,key=options.get),flush=True)


def finish_adaptation(args,d):
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        p=args.output/'extension_locked'/f'{part}.json'
        if not p.exists():complete_count_adaptation(args,d,part,outer,read(args.output/'tabm_locked'/f'{part}.json'))


def final(args,d):
    for part,_,_ in outer_partitions(d):assert (args.output/'extension_locked'/f'{part}.json').exists()
    for part,outer,test in outer_partitions(d):
        if args.parts and part not in args.parts:continue
        lock=read(args.output/'extension_locked'/f'{part}.json')
        for family in ['tabm','count_adapted']:
            choice=lock[family]
            if choice is None:continue
            for response in [True,False]:
                name=family+('_solvai' if response else '_structure');values=[]
                for seed in SEEDS:
                    path=args.output/'extension_final'/part/f'{name}_{seed}.json'
                    if path.exists():values.append(np.array(read(path)['predictions']));continue
                    r,cp,pred=train(d,choice['candidate'],outer,np.array([],dtype=int),seed,response,choice['refit_schedule'],part=='external')
                    y=pred(test);r.update(model=name,partition=part,test_ids=d.frame.id.iloc[test].tolist(),predictions=y.tolist())
                    if cp is not None:
                        cp_path=args.output/'extension_checkpoints'/f'{name}_{seed}.pt';cp_path.parent.mkdir(exist_ok=True)
                        torch.save(cp,cp_path);delta=float(np.max(np.abs(replay(d,cp,test)-y)));assert delta<1e-4
                        r.update(checkpoint_sha256=sha(cp_path),replay_max_error=delta)
                    dump(path,r);values.append(y);del cp,pred
                    print('EXT_REFIT',part,name,seed,flush=True)
                x=d.frame.iloc[test].copy();x['model']=name;x['prediction']=np.mean(values,axis=0)
                for seed,pred in zip(SEEDS,values):x[f'seed_{seed}']=pred
                p=args.output/'extension_predictions'/f'{name}_{part}.csv';p.parent.mkdir(exist_ok=True);x.to_csv(p,index=False)


def prepare(args,d):
    p=args.output/'extension_protocol.json';obj={'configs':configs(),'hashes':hashes(),'source_weights':d.source_hashes,
      'core_protocol_sha256':sha(args.output/'protocol.json'),'versions':{'tabm':'0.0.3','rtdl_num_embeddings':'0.0.12'},
      'before_any_final_test_evaluation':not (args.output/'predictions').exists(),'selection':'same nested equal-source metric; see EXTENSION.md'}
    assert obj['before_any_final_test_evaluation']
    if p.exists():assert read(p)==obj
    else:dump(p,obj)
    print('EXT_FROZEN',sha(p),flush=True)


def smoke(args,d):
    outer=np.flatnonzero(d.parts!='external');tr,va=split_inner(d,outer,0);reports=[]
    for k in [0,2,4,6,8,10]:
        c={**configs()[k],'cap':2,'schedule':'cosine'}
        r,cp,pred=train(d,c,tr[:128],va[:32],11,keep=True);a=pred(va[:32]);b=replay(d,cp,va[:32][::-1])[::-1]
        assert np.max(np.abs(a-b))<1e-4
        reports.append({'config':c['id'],'seconds':r['seconds'],'replay_and_permutation':True})
        del cp,pred
    delta=torch.tensor([[1.,-1.]])
    assert delta.square().mean().item()==1. and delta.mean(1).square().mean().item()==0.
    dump(args.output/'extension_smoke.json',{'tests':reports,'member_loss_test':True})
    print(reports,flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','smoke','search','adapt','final'])
    ap.add_argument('--release',type=Path,default=Path('RELEASE_ROOT'))
    ap.add_argument('--previous',type=Path,default=Path('BASE_INPUTS'))
    ap.add_argument('--output',type=Path,default=HERE/'results');ap.add_argument('--parts',nargs='+');args=ap.parse_args()
    torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    d=Dataset(args.release,args.previous)
    {'prepare':prepare,'smoke':smoke,'search':search,'adapt':finish_adaptation,'final':final}[args.stage](args,d)


if __name__=='__main__':main()
