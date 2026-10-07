"""Archived neural endpoint implementation; see ../README.md for data requirements."""
from __future__ import annotations
import hashlib
import json
import math
import random
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import Ridge

SEEDS = [11, 29, 47, 71, 101]
FAMILIES = ['resnet', 'dual', 'ridge_residual', 'embedding', 'graph']
INPUT_HASHES = {
    'molecules.csv': 'ebd1787cabe2326193d795a4fa478cfaf7e1df0b3c2beb2c533c0bb589153c20',
    'features.npz': '90ac6a5f0f55ba2aabe41a2ac4c82d98721777d0b985f8023d3d2822fc747741',
    'molformer.npy': 'b47560bb6da54f802af295231692253180defb14f149ab3cce57479baf1acf7e',
    'unimol.npy': 'd4b62d478036e4659dd76188e067289f3641d3ae51f9a2e0c4be4ae422f386d0',
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(p, obj):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_suffix(p.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False, default=str) + '\n')
    temp.replace(p)


def candidates():
    # Balanced, fixed alternatives, not outcomes from the external cohort.
    settings = [
        dict(width=32, depth=1, dropout=0., decay=.001, lr=.003, loss='mse', schedule='plateau'),
        dict(width=64, depth=2, dropout=.1, decay=.01, lr=.001, loss='huber', schedule='plateau'),
        dict(width=128, depth=2, dropout=.1, decay=.001, lr=.0003, loss='mse', schedule='cosine'),
        dict(width=256, depth=3, dropout=.2, decay=.01, lr=.001, loss='huber', schedule='cosine'),
    ]
    result = []
    for family in FAMILIES[:3]:
        for representation in ['binary_log', 'count_log', 'count_drop']:
            for setting in settings:
                result.append(dict(family=family, representation=representation, **setting))
    for encoder in ['molformer', 'unimol']:
        for setting in settings:
            result.append(dict(family='embedding', representation='count_log', encoder=encoder, **setting))
    for readout in ['sum_mean', 'atomic_global']:
        for j, setting in enumerate(settings):
            # Real Chemprop bond message passing; not an MLP labelled as a GNN.
            result.append(dict(family='graph', representation='count_log', readout=readout,
                               message_dim=[64, 128, 128, 256][j], message_depth=[3, 3, 5, 4][j],
                               **setting))
    for i, c in enumerate(result):
        c.update(id=f'c{i:02d}', cap=600, batch=128)
    assert len(result) == 52
    return result


def seed_all(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)


class Dataset:
    def __init__(self, release, previous):
        self.release, self.previous = Path(release), Path(previous)
        b = self.release / 'results/journal_20261002'
        for n, h in INPUT_HASHES.items(): assert sha(b / n) == h, n
        self.frame = pd.read_csv(b/'molecules.csv', float_precision='round_trip')
        assert len(self.frame) == 1585 and self.frame.id.is_unique
        assert self.frame.part.value_counts().to_dict() == {'public':1280, 'arrow':85, 'external':220}
        p = json.loads((self.previous/'results/protocol.json').read_text())
        self.columns = p['feature_columns'][:2265]
        assert all(c.startswith('morgan2__') for c in self.columns[217:])
        self.source_hashes = p['source_weights']
        for n, h in self.source_hashes.items(): assert sha(self.release/'models/final'/n) == h, n
        z = np.load(b/'features.npz')
        self.struct, self.response = z['structure'].astype(float), z['response'].astype(float)
        self.embeddings = {k:np.load(b/f'{k}.npy') for k in ['molformer','unimol']}
        self.smiles = self.frame.smiles.tolist()
        self.mols = [Chem.MolFromSmiles(s) for s in self.smiles]
        assert all(m is not None for m in self.mols)
        self.size = np.array([m.GetNumHeavyAtoms() for m in self.mols])
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        self.counts = np.array([gen.GetCountFingerprintAsNumPy(m) for m in self.mols], dtype=float)
        self.graphs = None
        self.y = self.frame.y.to_numpy(float)
        self.weights = self.frame.weight.to_numpy(float)
        self.parts = self.frame.part.to_numpy()
        self.stress_data = None

    def arrays(self, c, response=True):
        rd = self.struct[:,:217].copy()
        ipc = self.columns.index('rdkit__Ipc')
        if c['representation'].endswith('drop'): rd[:,ipc] = 0.
        else: rd[:,ipc] = np.log1p(np.maximum(rd[:,ipc], 0.))
        # Bertz and Kappa descriptors can span large ranges. These fixed transforms
        # leave physical counts, molecular weight and surface areas uncompressed.
        for name in ['BertzCT','Kappa1','Kappa2','Kappa3']:
            j = self.columns.index('rdkit__'+name)
            rd[:,j] = np.sign(rd[:,j])*np.log1p(np.abs(rd[:,j]))
        if c['family'] == 'graph': x = rd
        else:
            fp = self.struct[:,217:] if c['representation'].startswith('binary') else self.counts
            x = np.column_stack([rd,fp])
        if c['family'] == 'embedding': x = np.column_stack([x,self.embeddings[c['encoder']]])
        n_structure = x.shape[1]
        if response: x = np.column_stack([x,self.response])
        binary = np.zeros(x.shape[1], dtype=bool)
        if c['representation'].startswith('binary') and c['family'] != 'graph': binary[217:2265] = True
        return x, n_structure, binary

    def get_graphs(self):
        if self.graphs is None:
            from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
            featurizer = SimpleMoleculeMolGraphFeaturizer()
            self.graphs = [featurizer(m) for m in self.mols]
        return self.graphs


def outer_partitions(d):
    f = d.frame
    for i in range(5):
        yield f'arrow_{i}', np.flatnonzero((f.part=='public')|((f.part=='arrow')&(f.fold!=i))), np.flatnonzero((f.part=='arrow')&(f.fold==i))
    yield 'external', np.flatnonzero(f.part!='external'), np.flatnonzero(f.part=='external')


def split_inner(d, outer, index):
    rng = np.random.default_rng(640601 + index)
    val = []
    for source in ['public','arrow']:
        eligible = outer[d.parts[outer] == source]
        n = max(1, round(.15*len(eligible)))
        if index == 3:
            order = sorted(eligible, key=lambda i:(-int(d.size[i]),str(d.frame.id.iloc[i])))
        else: order = rng.permutation(eligible)
        val.extend(order[:n])
    val = np.array(sorted(val),dtype=int); tr = np.setdiff1d(outer,val)
    assert len(np.intersect1d(tr,val)) == 0 and set(tr)|set(val) == set(outer)
    assert not np.any(d.parts[outer]=='external')
    return tr,val


class Preprocessor:
    def fit(self, x, boundary, binary):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            self.median = np.nan_to_num(np.nanmedian(np.where(np.isfinite(x),x,np.nan),axis=0))
        x = np.where(np.isfinite(x),x,self.median)
        # Keep all response dimensions even if one happens to be constant.
        self.keep = (np.ptp(x,axis=0)>1e-12) | (np.arange(x.shape[1])>=boundary)
        self.boundary = int(self.keep[:boundary].sum())
        self.ridge_indices = np.flatnonzero((np.flatnonzero(self.keep)<217)|(np.flatnonzero(self.keep)>=boundary))
        a = x[:,self.keep]
        self.mean, self.scale = a.mean(0), a.std(0)
        self.scale[self.scale<1e-8] = 1.
        b = binary[self.keep]; self.mean[b]=0.; self.scale[b]=1.
        return self

    def transform(self,x):
        a = np.where(np.isfinite(x),x,self.median)[:,self.keep]
        a = ((a-self.mean)/self.scale).astype(np.float32)
        assert np.isfinite(a).all()
        return a

    def state(self): return dict(vars(self))


class Block(nn.Module):
    def __init__(self,w,p):
        super().__init__()
        self.layers=nn.Sequential(nn.LayerNorm(w),nn.Linear(w,w*2),nn.GELU(),nn.Dropout(p),nn.Linear(w*2,w),nn.Dropout(p))
    def forward(self,x): return x+self.layers(x)


def tower(inp,w,depth,drop):
    return nn.Sequential(nn.Linear(inp,w),nn.GELU(),*[Block(w,drop) for _ in range(depth)])


class Head(nn.Module):
    def __init__(self,config,nin,boundary):
        super().__init__(); self.c=config; self.boundary=boundary
        w,p,depth=config['width'],config['dropout'],config['depth']
        self.has_response = nin>boundary
        self.graph = config['family']=='graph'
        if self.graph:
            from chemprop.nn import BondMessagePassing
            h=config['message_dim']
            self.mp=BondMessagePassing(d_h=h, depth=config['message_depth'], dropout=p)
            if config['readout']=='atomic_global':
                self.atomic=nn.Sequential(nn.Linear(h,w),nn.GELU(),nn.Linear(w,1))
                graphdim=h+1
            else: graphdim=h*2+1
            self.struct=tower(boundary+graphdim,w,depth,p)
        elif config['family'] in ['dual','embedding']:
            self.struct=tower(boundary,w,depth,p)
        else:
            self.struct=tower(nin,w,depth,p)
        self.separate = self.graph or config['family'] in ['dual','embedding']
        if self.separate and self.has_response:
            self.response=tower(nin-boundary,32,1,p)
            outdim=w+32
        else: outdim=w
        self.output=nn.Sequential(nn.Linear(outdim,w),nn.GELU(),nn.Linear(w,1))
        nn.init.zeros_(self.output[-1].bias)
        nn.init.normal_(self.output[-1].weight,std=.01)
        if self.graph and config['readout']=='atomic_global':
            nn.init.zeros_(self.atomic[-1].bias); nn.init.normal_(self.atomic[-1].weight,std=.01)

    def forward(self,x,g=None):
        atom_out=0.
        if self.graph:
            h=self.mp(g); n=x.shape[0]
            sums=torch.zeros(n,h.shape[1],device=x.device,dtype=h.dtype).index_add_(0,g.batch,h)
            count=torch.bincount(g.batch,minlength=n).to(h.dtype).unsqueeze(1)
            means=sums/count
            if self.c['readout']=='atomic_global':
                a=self.atomic(h)
                atom_out=torch.zeros(n,1,device=x.device,dtype=h.dtype).index_add_(0,g.batch,a).squeeze(1)/10.
                info=torch.cat([means,torch.log1p(count)],1)
            else: info=torch.cat([sums/10.,means,torch.log1p(count)],1)
            z=self.struct(torch.cat([x[:,:self.boundary],info],1))
        else: z=self.struct(x[:,:self.boundary] if self.separate else x)
        if self.separate and self.has_response: z=torch.cat([z,self.response(x[:,self.boundary:])],1)
        return self.output(z).squeeze(1)+atom_out


def scores(d,idx,pred):
    err=np.abs(pred-d.y[idx])
    per={p:float(err[d.parts[idx]==p].mean()) for p in ['public','arrow'] if np.any(d.parts[idx]==p)}
    return {'score':float(np.mean(list(per.values()))),'pooled_mae':float(err.mean()),'by_source':per}


def fit(d,c,train,valid,seed,with_response=True,schedule_override=None,checkpoint=False):
    seed_all(seed); start=time.monotonic()
    assert not np.any(d.parts[train]=='external')
    assert not np.any(d.parts[valid]=='external')
    assert not np.intersect1d(train,valid).size
    x,b,binary=d.arrays(c,with_response)
    prep=Preprocessor().fit(x[train],b,binary); xx=prep.transform(x)
    mean=float(d.y[train].mean()); scale=float(d.y[train].std())
    yy=(d.y-mean)/scale
    base=np.zeros(len(yy),dtype=float); ridge_state=None
    if c['family']=='ridge_residual':
        cols=prep.ridge_indices
        ridge=Ridge(alpha=10.,solver='cholesky').fit(xx[train][:,cols],yy[train],sample_weight=d.weights[train])
        base=ridge.predict(xx[:,cols]).astype(float)
        ridge_state={'columns':cols,'coef':ridge.coef_,'intercept':float(ridge.intercept_)}
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    xt=torch.as_tensor(xx,device=device); yt=torch.as_tensor((yy-base).astype(np.float32),device=device)
    wt=torch.as_tensor((d.weights/d.weights[train].mean()).astype(np.float32),device=device)
    model=Head(c,xx.shape[1],prep.boundary).to(device)
    optimizer=torch.optim.AdamW([
        {'params':[p for p in model.parameters() if p.ndim>=2], 'weight_decay':c['decay']},
        {'params':[p for p in model.parameters() if p.ndim<2], 'weight_decay':0.}],lr=c['lr'])
    plateau=torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,factor=.5,patience=25,threshold=0.,min_lr=1e-6)
    graphs=d.get_graphs() if model.graph else None
    graph_cache={}
    def batch_graph(ids):
        if graphs is None: return None
        key=tuple(ids)
        if key not in graph_cache:
            from chemprop.data import BatchMolGraph
            g=BatchMolGraph([graphs[i] for i in ids]); g.to(device); graph_cache[key]=g
        return graph_cache[key]
    # Graph minibatches are fixed per fit and shuffled in order each epoch, so cached
    # graph memory is bounded. No cross-molecule operation mixes representations.
    rng=np.random.default_rng(seed)
    train_batches=np.array_split(rng.permutation(train),math.ceil(len(train)/c['batch']))
    def predict_ids(ids):
        model.eval(); preds=[]
        with torch.inference_mode():
            for batch in np.array_split(ids,max(1,math.ceil(len(ids)/256))):
                if not len(batch): continue
                preds.extend(model(xt[batch],batch_graph(batch)).cpu().numpy().astype(float))
        return (np.asarray(preds)+base[ids])*scale+mean
    best=float('inf'); best_epoch=0; best_pred=None; trace=[]
    cap=len(schedule_override) if schedule_override is not None else c['cap']
    for epoch in range(1,cap+1):
        if schedule_override is not None: lr=float(schedule_override[epoch-1])
        elif c['schedule']=='cosine':
            lr=c['lr']*min(1.,epoch/10.)*(.01+.99*(1+math.cos(math.pi*(epoch-1)/cap))/2)
        else: lr=optimizer.param_groups[0]['lr']
        for group in optimizer.param_groups: group['lr']=lr
        model.train()
        if model.graph:
            batches=[train_batches[j] for j in rng.permutation(len(train_batches))]
        else: batches=np.array_split(rng.permutation(train),math.ceil(len(train)/c['batch']))
        for ids in batches:
            optimizer.zero_grad(set_to_none=True)
            delta=model(xt[ids],batch_graph(ids))-yt[ids]
            if c['loss']=='mse': losses=delta.square()
            else: losses=nn.functional.huber_loss(delta,torch.zeros_like(delta),delta=1.,reduction='none')
            loss=(losses*wt[ids]).mean()
            if not torch.isfinite(loss): raise FloatingPointError(f'Nonfinite loss {c["id"]} {epoch}')
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5.); optimizer.step()
        if len(valid):
            pred=predict_ids(valid); sc=scores(d,valid,pred)
            row={'epoch':epoch,'lr':lr,**sc}
            if epoch==1 or epoch%10==0: row['train']=scores(d,train,predict_ids(train))
            trace.append(row)
            if epoch%100==0:
                print('EPOCH',c['id'],seed,epoch,'validation',round(sc['score'],5),'lr',round(lr,8),flush=True)
            if sc['score']<best:
                best=sc['score']; best_epoch=epoch; best_pred=pred.copy()
            if c['schedule']=='plateau': plateau.step(sc['score'])
            # Full cosine schedules are deliberately completed, not cut before annealing.
            if c['schedule']=='plateau' and epoch>=150 and epoch-best_epoch>=100:
                break
    record={'config':c,'seed':seed,'with_response':with_response,'epochs_run':epoch,
            'best_epoch':best_epoch if len(valid) else epoch,'best_score':best if len(valid) else None,
            'seconds':time.monotonic()-start,'trace':trace,'parameters':sum(p.numel() for p in model.parameters()),
            'train_ids':d.frame.id.iloc[train].tolist(),'validation_ids':d.frame.id.iloc[valid].tolist(),
            'validation_predictions':best_pred.tolist() if best_pred is not None else [],
            'train_final':scores(d,train,predict_ids(train))}
    cp=None
    if checkpoint:
        state={k:v.detach().cpu() for k,v in model.state_dict().items()}
        cp={'config':c,'state':state,'prep':prep.state(),'mean':mean,'scale':scale,'ridge':ridge_state,
            'nin':xx.shape[1],'boundary':prep.boundary,'with_response':with_response,'seed':seed}
    return record,cp,predict_ids


def schedule_from_records(records):
    epochs=max(1,round(float(np.median([r['best_epoch'] for r in records]))))
    schedules=[]
    for r in records:
        lr=[t['lr'] for t in r['trace']]
        schedules.append([lr[min(i,len(lr)-1)] for i in range(epochs)])
    return np.median(schedules,axis=0).tolist()


def checkpoint_predict(d,cp,ids):
    x,_,_=d.arrays(cp['config'],cp['with_response'])
    prep=Preprocessor(); prep.__dict__.update(cp['prep']); xx=prep.transform(x)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model=Head(cp['config'],cp['nin'],cp['boundary']).to(device)
    model.load_state_dict(cp['state']); model.eval(); values=[]
    with torch.inference_mode():
        for b in np.array_split(ids,max(1,math.ceil(len(ids)/128))):
            g=None
            if model.graph:
                from chemprop.data import BatchMolGraph
                gs=d.get_graphs(); g=BatchMolGraph([gs[i] for i in b]);g.to(device)
            z=model(torch.as_tensor(xx[b],device=device),g).cpu().numpy().astype(float)
            if cp['ridge'] is not None:
                r=cp['ridge']; z+=xx[b][:,r['columns']]@r['coef']+r['intercept']
            values.extend(z*cp['scale']+cp['mean'])
    return np.asarray(values)
