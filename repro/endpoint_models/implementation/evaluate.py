"""Post-selection verification and evaluation; never selects or refits a neural model."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from core import Dataset, FAMILIES, SEEDS, checkpoint_predict, dump, outer_partitions, sha, split_inner


def prepare_stress(args,d):
    b=args.previous/'results'
    grid=pd.read_csv(b/'stress_grid.csv',float_precision='round_trip')
    assert len(grid)==56
    grid['canonical_smiles']=[Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in grid.smiles]
    grid.to_csv(args.output/'stress_grid.csv',index=False)
    sys=__import__('sys');sys.path.insert(0,str(args.release))
    from solv_ai.features import descriptor_frame
    dest=args.output/'stress_features.npz'
    if not dest.exists():
        struct=descriptor_frame(grid.smiles.tolist())[d.columns].to_numpy(np.float32)
        response=np.load(b/'stress_features.npz')['response']
        assert response.shape==(56,15)
        mols=[Chem.MolFromSmiles(s) for s in grid.smiles]
        generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        counts=np.array([generator.GetCountFingerprintAsNumPy(m) for m in mols],dtype=np.float32)
        assert np.array_equal(counts>0,struct[:,217:]>0)
        np.savez_compressed(dest,structure=struct,response=response,counts=counts)
    # Reuse known molecule encodings from the frozen cache. Generate only new identities.
    existing={Chem.MolToSmiles(m):i for i,m in enumerate(d.mols)}
    metadata={}
    for name,width in [('molformer',768),('unimol',512)]:
        path=args.output/f'stress_{name}.npy'
        x=np.full((56,width),np.nan,dtype=np.float32)
        unseen=[]
        for i,s in enumerate(grid.canonical_smiles):
            if s in existing: x[i]=d.embeddings[name][existing[s]]
            else: unseen.append(i)
        smiles=grid.canonical_smiles.iloc[unseen].tolist()
        if path.exists():
            a=np.load(path);assert a.shape==(56,width) and np.isfinite(a).all()
            metadata[name]={'cached':56-len(unseen),'new':len(unseen),'sha256':sha(path),'shape':list(a.shape),
                            'labels_used':False,'canonical_smiles_used':True,'batch_size':8}
            continue
        print('STRESS_ENCODER',name,'cached',56-len(unseen),'new',len(unseen),flush=True)
        if name=='molformer':
            from transformers import AutoModel,AutoTokenizer
            model_id='ibm-research/MoLFormer-XL-both-10pct'
            rev='7b12d946c181a37f6012b9dc3b002275de070314'
            tokenizer=AutoTokenizer.from_pretrained(model_id,revision=rev,trust_remote_code=True,local_files_only=True)
            encoder=AutoModel.from_pretrained(model_id,revision=rev,trust_remote_code=True,local_files_only=True).eval().cuda()
            token_lengths=[len(tokenizer(s)['input_ids']) for s in smiles]
            assert max(token_lengths,default=0)<=256,'Unexpected SMILES truncation'
            with torch.inference_mode():
                for start in range(0,len(smiles),8):
                    batch=tokenizer(smiles[start:start+8],padding=True,truncation=True,max_length=256,return_tensors='pt')
                    batch={k:v.cuda() for k,v in batch.items()}
                    h=encoder(**batch).last_hidden_state;mask=batch['attention_mask'].unsqueeze(-1).to(h.dtype)
                    x[np.array(unseen[start:start+8])]=((h*mask).sum(1)/mask.sum(1).clamp_min(1)).cpu().numpy()
            del encoder
        else:
            from unimol_tools import UniMolRepr
            encoder=UniMolRepr(data_type='molecule',batch_size=8,remove_hs=False,
                               model_name='unimolv1',model_size='84m',use_cuda=True)
            x[unseen]=np.asarray(encoder.get_repr(smiles))
            del encoder
        assert np.isfinite(x).all()
        np.save(path,x)
        metadata[name]={'cached':56-len(unseen),'new':len(unseen),'sha256':sha(path),'shape':list(x.shape),
                        'labels_used':False,'canonical_smiles_used':True,'batch_size':8}
        torch.cuda.empty_cache()
    dump(args.output/'stress_feature_metadata.json',{'source_weights':d.source_hashes,'features_sha256':sha(dest),
                                                    'encoders':metadata,'previous_response_sha256':sha(b/'stress_features.npz')})


def load_stress(args,d):
    z=np.load(args.output/'stress_features.npz');grid=pd.read_csv(args.output/'stress_grid.csv')
    s=copy.copy(d)
    s.frame=grid;s.smiles=grid.smiles.tolist();s.mols=[Chem.MolFromSmiles(x) for x in s.smiles]
    s.struct=z['structure'].astype(float);s.response=z['response'].astype(float);s.counts=z['counts'].astype(float)
    s.embeddings={n:np.load(args.output/f'stress_{n}.npy') for n in ['molformer','unimol']}
    if type(d).__name__=='CorrectedDataset':
        s.embeddings['molformer']=np.load(args.output/'stress_molformer_deterministic.npy')
    s.graphs=None
    return s,grid


def stress(args,d):
    s,grid=load_stress(args,d);rows=[]
    for family in FAMILIES:
        for condition in ['solvai','structure']:
            model=f'{family}_{condition}';values=[]
            for seed in SEEDS:
                cp=torch.load(args.output/'checkpoints'/f'{model}_{seed}.pt',map_location='cpu',weights_only=False)
                values.append(checkpoint_predict(s,cp,np.arange(len(grid))))
            out=grid.copy();out['model']=model;out['prediction']=np.mean(values,axis=0)
            out['seed_std']=np.std(values,axis=0,ddof=1)
            for seed,v in zip(SEEDS,values):out[f'seed_{seed}']=v
            assert np.isfinite(out.prediction).all()
            for n in range(13):
                if n==0:
                    same=out[(out.n==0)&out.series.isin(['glycine','alanine'])].prediction.to_numpy()
                    assert abs(same[0]-same[1])<1e-4
            rows.append(out)
            print('STRESS_DONE',model,float(out.prediction.min()),float(out.prediction.max()),flush=True)
    pd.concat(rows,ignore_index=True).to_csv(args.output/'stress_predictions.csv',index=False)


def paired(a,b):
    a=a.sort_values('id');b=b.sort_values('id')
    assert a.id.tolist()==b.id.tolist() and np.array_equal(a.y,b.y)
    delta=(a.prediction-a.y).abs().to_numpy()-(b.prediction-b.y).abs().to_numpy()
    rng=np.random.default_rng(20260828)
    boot=np.concatenate([delta[rng.integers(0,len(delta),(1000,len(delta)))].mean(axis=1) for _ in range(100)])
    lo,hi=np.quantile(boot,[.025,.975])
    return float(delta.mean()),float(lo),float(hi)


def summarize(args,d):
    rows=[]
    for part,_,test in outer_partitions(d):
        lock=json.loads((args.output/'locked'/f'{part}.json').read_text())
        for family in FAMILIES:
            for condition in ['solvai','structure']:
                model=f'{family}_{condition}'
                x=pd.read_csv(args.output/'predictions'/f'{model}_{part}.csv',float_precision='round_trip')
                assert x.id.tolist()==d.frame.id.iloc[test].tolist()
                assert np.max(np.abs(x.prediction-x[[f'seed_{s}' for s in SEEDS]].mean(axis=1)))<1e-12
                rows.append(x)
                if family==lock['overall_family']:
                    x=x.copy();x['model']='selected_'+condition;rows.append(x)
    previous=pd.read_csv(args.previous/'results/all_test_predictions.csv',float_precision='round_trip')
    previous.loc[~previous.model.str.startswith('tree_'),'model']='pilot_'+previous.loc[~previous.model.str.startswith('tree_'),'model']
    rows.append(previous)
    f=pd.concat(rows,ignore_index=True)
    assert not f.duplicated(['id','model']).any()
    assert f.groupby('model').size().eq(305).all()
    assert np.isfinite(f.prediction).all()
    f.to_csv(args.output/'all_test_predictions.csv',index=False)
    metrics=[];groups={}
    for model,g in f.groupby('model'):
        for cohort,mask in [('ARROW-85',g.part=='arrow'),('External-220',g.part=='external'),('Strict-97',(g.part=='external')&g.strict)]:
            a=g[mask];err=(a.prediction-a.y).to_numpy();groups[cohort,model]=a
            available=[f'seed_{s}' for s in SEEDS if f'seed_{s}' in a and a[f'seed_{s}'].notna().all()]
            seed_maes=[float(np.abs(a[s]-a.y).mean()) for s in available]
            metrics.append({'cohort':cohort,'model':model,'n':len(a),'mae':float(np.abs(err).mean()),
                            'rmse':float(np.sqrt(np.mean(err**2))),'median_ae':float(np.median(np.abs(err))),
                            'seed_mae_std':float(np.std(seed_maes,ddof=1)),'seed_count':len(available)})
    pd.DataFrame(metrics).to_csv(args.output/'metrics.csv',index=False)
    comparisons=[]
    for cohort in ['ARROW-85','External-220','Strict-97']:
        for family in FAMILIES+['selected']:
            for ref in ['tree_solvai',family+'_structure','pilot_resnet_solvai','pilot_linear_resnet_solvai']:
                name=family+'_solvai';delta,lo,hi=paired(groups[cohort,name],groups[cohort,ref])
                comparisons.append({'cohort':cohort,'candidate':name,'reference':ref,'delta':delta,'low':lo,'high':hi})
    pd.DataFrame(comparisons).to_csv(args.output/'paired_comparisons.csv',index=False)
    print(pd.DataFrame(metrics).pivot(index='model',columns='cohort',values='mae').to_string(),flush=True)


def audit(args,d):
    protocol=json.loads((args.output/'protocol.json').read_text())
    for n,h in protocol['source_hashes'].items(): assert sha(Path(__file__).parent/n)==h,n
    records=0;cap_hits=[];choice_rows=[]
    for part,outer,test in outer_partitions(d):
        outer_ids=set(d.frame.id.iloc[outer]); test_ids=set(d.frame.id.iloc[test])
        mapping={s:i for i,s in enumerate(d.frame.id)}
        for path in (args.output/'inner'/part).glob('*.json'):
            r=json.loads(path.read_text());assert not r.get('failed'),path
            tr=set(r['train_ids']);va=set(r['validation_ids'])
            assert tr|va==outer_ids and not tr&va and not (tr|va)&test_ids
            assert r['source_hashes']==protocol['source_hashes']
            trace=r['trace'];best=min(trace,key=lambda t:t['score'])
            assert best['epoch']==r['best_epoch'] and abs(best['score']-r['best_score'])<1e-12
            idx=np.array([mapping[i] for i in r['validation_ids']]); pred=np.array(r['validation_predictions'])
            scores=[np.abs(pred[d.parts[idx]==p]-d.y[idx][d.parts[idx]==p]).mean() for p in ['public','arrow']]
            assert abs(np.mean(scores)-r['best_score'])<1e-10
            assert len(trace)==r['epochs_run']
            if r['best_epoch']>=.95*r['config']['cap']:cap_hits.append({'path':path.name,'part':part,'best_epoch':r['best_epoch']})
            records+=1
        lock=json.loads((args.output/'locked'/f'{part}.json').read_text())
        assert lock['source_hashes']==protocol['source_hashes']
        assert lock['overall_family']==min(lock['winners'],key=lambda a:(a['mean_score'],a['family']))['family']
        for w in lock['winners']:
            choice_rows.append({'partition':part,**{k:w[k] for k in ['family','use_long','mean_score','chosen_epochs']},
                                'config':w['candidate'],'refit_epochs':len(w['refit_schedule']),
                                'long_minus_normal':w['long_minus_normal']})
        for path in (args.output/'final_runs'/part).glob('*.json'):
            r=json.loads(path.read_text())
            assert set(r['train_ids'])==outer_ids and not r['validation_ids']
            assert set(r['test_ids'])==test_ids and r['source_hashes']==protocol['source_hashes']
            if part=='external':assert sha(args.output/'checkpoints'/f'{r["model"]}_{r["seed"]}.pt')==r['checkpoint_sha256']
    dump(args.output/'audit.json',{'status':'pass','inner_fits':records,'cap_near_best_records':cap_hits,
                                  'choices':choice_rows,'source_and_target_disjointness':True,
                                  'validation_scores_recomputed_from_labels':True,
                                  'no_physical_accuracy_claim_for_unlabelled_homologues':True})
    print('AUDIT_PASS',records,'inner fits;',len(cap_hits),'near-cap best epochs',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare_stress','stress','summarize','audit'])
    ap.add_argument('--release',type=Path,default=Path('RELEASE_ROOT'))
    ap.add_argument('--previous',type=Path,default=Path('BASE_INPUTS'))
    ap.add_argument('--output',type=Path,default=Path(__file__).parent/'results')
    args=ap.parse_args();torch.set_num_threads(2)
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    d=Dataset(args.release,args.previous)
    globals()[args.stage](args,d)


if __name__=='__main__':main()
