#!/usr/bin/env python3
"""Build the journal figure suite locally from frozen data and editable vectors."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import math
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from journal_vector_art import Art, ROOT, ASSETS, INK, MID, GRID, BLUE, TEAL, AMBER, PURPLE, ROSE, COLORS, COUNTS

OUT = ROOT / 'paper/figures/journal'
SI = ROOT / 'paper/supplementary/figures'
JR = ROOT / 'results/journal_20261002'
ER = ROOT / 'results/endpoint_models_20261006'
COHORTS = ['ARROW-85', 'External-220', 'Strict-97']

plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,
    'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':.6,
    'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':MID,
    'axes.labelcolor':INK,'xtick.color':MID,'ytick.color':MID,'text.color':INK,
    'svg.fonttype':'none','pdf.fonttype':42,'svg.hashsalt':'solvai-local-journal'})

from journal_style import apply as apply_style
apply_style(plt)


def load_data():
    endpoint=pd.read_parquet(ROOT/'results/confirmatory/standardized_exclusion_endpoint_predictions.parquet')
    return dict(metrics=json.loads((ROOT/'results/paper_metrics.json').read_text()),
        primary=endpoint.loc[endpoint.partition.eq('standardized_exclusion_primary')],
        repeats=endpoint.loc[endpoint.partition.eq('standardized_exclusion_repeat')],
        zero=endpoint.loc[endpoint.partition.eq('standardized_exclusion_zero_arrow')],
        oldpairs=pd.read_csv(ROOT/'results/confirmatory/confirmatory_paired_comparisons.csv'),
        separation=pd.read_csv(ROOT/'results/confirmatory/standardized_exclusion_global_separation_metrics.csv'),
        modern=pd.read_csv(ER/'corrected_original_metrics.csv'), pairs=pd.read_csv(ER/'corrected_original_comparisons.csv'),
        strata=pd.read_csv(JR/'similarity_strata.csv'))


def mae(data,model,cohort):
    f=data['modern']
    return float(f.loc[f.model.eq(model)&f.cohort.eq(cohort),'mae'].item())


def base_plot(w,h,rect=(.18,.20,.78,.74)):
    fig=plt.figure(figsize=(w*.25/25.4,h*.25/25.4))
    ax=fig.add_axes(rect)
    return fig,ax


def grid(ax,axis='x'):
    ax.set_axisbelow(True);ax.grid(axis=axis,color=GRID,lw=.5)
    ax.tick_params(width=.6,length=3)


def source_glyph(c,i,x,y,w=66,h=55):
    if i in (0,3,4):
        c.shell(x,y,w,h,explicit=False)
    elif i==2:
        c.shell(x,y,w,h,explicit=True)
    elif i==1:
        cx,cy=x+w/2,y+h/2
        c.molecule(cx-16,cy-14,32,29,hydrogens=False)
        for j,label in enumerate(['E','S','A','B','L']):
            angle=-math.pi/2+j*2*math.pi/5
            xx,yy=cx+34*math.cos(angle),cy+27*math.sin(angle)
            c.line(cx,cy,xx,yy,AMBER,.7)
            c.circle(xx,yy,6.5,'white')
            c.text(xx,yy+3.5,label,10,AMBER,600,'middle')
    else:
        for j,m in enumerate(c.conformers):
            c.molecule(x+j*23,y+12,26,33,mol=m,hydrogens=False,angle=-20+j*10)


def overview(data,preview):
    from journal_teaser import overview as draw
    return draw(data,preview)


def learning(data,preview):
    from journal_teaser import learning as draw
    return draw(data,preview)


def matched(data,preview):
    c=Art(720,'Matched tests of molecule-aligned solvent-response information')
    c.title(10,25,'a','Change only the response information')
    c.text(34,47,'Same structures, labels, folds, weights, learner and seeds',12,MID)
    shuffle=data['oldpairs'].loc[data['oldpairs'].analysis.eq('aligned_vs_mean_shuffle_primary_-1'),'reference_mae'].item()
    conditions=[('Structure only','none',mae(data,'structure','ARROW-85'),MID),
        ('Aligned responses','aligned',mae(data,'solvai','ARROW-85'),TEAL),
        ('Shuffled responses','shuffled',shuffle,AMBER)]
    for i,(name,mode,value,color) in enumerate(conditions):
        x=17+i*238
        c.text(x+106,79,name,13,color,600,'middle')
        for j in range(3):
            yy=102+j*16;c.circle(x+9,yy,3,INK)
            c.feature_stripes(x+20,yy-4,25,8)
            if mode!='none':
                k=([2,0,1][j] if mode=='shuffled' else j)
                c.line(x+48,yy,x+76,102+k*16,color,.8)
                for t in range(5):c.rect(x+79+t*6,98+k*16,4.5,8,COLORS[t],opacity=[.4,.7,1][k])
        c.arrow(x+117,117,x+141,117,color)
        c.text(x+174,125,f'{value:.3f}',24,color,600,'middle')
        c.text(x+174,144,'MAE',10.5,MID,400,'middle')
    c.line(10,163,710,163,GRID)
    c.title(10,192,'b','Build the response set')
    c.title(389,192,'c','Resolve molecule-level changes')
    fig,ax=base_plot(342,225,(.17,.23,.77,.69))
    methods=data['metrics']['methods'];keys=['matched_structure_only','narrow_response','narrow_plus_smd','full_solvai']
    values=[methods[k]['mae_kcal_mol'] for k in keys]
    ax.plot(range(4),values,color=BLUE,lw=1.3)
    for i,v in enumerate(values):
        ax.scatter(i,v,s=35,color=MID if i==0 else TEAL if i==3 else BLUE,zorder=3)
        ax.text(i,v+.008,f'{v:.3f}',ha='center',fontsize=8.2,color=INK)
    pimd=methods['arrow_pimd8']['mae_kcal_mol'];ax.axhline(pimd,color=ROSE,lw=.9,ls='--')
    ax.text(.04,pimd-.012,'PIMD8',transform=ax.get_yaxis_transform(),
            ha='left',va='top',color=ROSE,fontsize=8)
    ax.set(xlim=(-.25,3.25),ylim=(.177,.333),xticks=range(4),
        xticklabels=['Structure','+ 8','+ SMD','+ ConfSolv'],ylabel='OOF MAE (kcal mol$^{-1}$)')
    grid(ax,'y');c.embed_plot(fig,10,210,342,225,'progression');plt.close(fig)
    fig,ax=base_plot(325,225,(.21,.23,.72,.69))
    p=data['primary'];a=p.loc[p.method.eq('A_structure_only')].sort_values('molecule_id');b=p.loc[p.method.eq('F_full_solvai')].sort_values('molecule_id')
    assert list(a.molecule_id)==list(b.molecule_id)
    aa,bb=a.absolute_error.to_numpy(),b.absolute_error.to_numpy();limit=max(aa.max(),bb.max())*1.05
    ax.fill_between([0,limit],[0,limit],color=TEAL,alpha=.065)
    ax.plot([0,limit],[0,limit],ls='--',color=MID,lw=.7)
    ax.scatter(aa,bb,s=17,color=np.where(bb<aa,TEAL,MID),edgecolors='white',lw=.35,zorder=3)
    ax.text(.97,.06,f'{int((bb<aa).sum())}/85 improved',transform=ax.transAxes,ha='right',fontsize=8,color=TEAL)
    ax.set(xlim=(0,limit),ylim=(0,limit),xlabel='Structure-only error (kcal mol$^{-1}$)',ylabel='SolvAI error (kcal mol$^{-1}$)')
    ax.set_aspect('equal',adjustable='box')
    c.embed_plot(fig,389,210,325,225,'molecules');plt.close(fig)
    c.title(10,477,'d','Which source blocks contribute?')
    c.title(389,477,'e','Repeat the complete partition')
    fig,ax=base_plot(342,222,(.34,.23,.60,.67));pairs=data['oldpairs'].set_index('analysis')
    keys=['primary_B_empirical_residual','primary_C_computation_core','primary_D_smd_water','primary_E_confsolv','primary_F_full_solvai']
    for i,key in enumerate(keys):
        r=pairs.loc[key];color=TEAL if i==4 else BLUE
        ax.errorbar(r.difference,i,xerr=[[r.difference-r.ci_low_95],[r.ci_high_95-r.difference]],fmt='o',ms=4.5,color=color,lw=1.2,capsize=2.5)
    ax.axvline(0,color=MID,lw=.8);ax.set(yticks=range(5),yticklabels=['Empirical +\ncorrected','Computed core','SMD','ConfSolv','Full 15'],ylim=(4.6,-.6),xlabel='ΔMAE vs structure (kcal mol$^{-1}$)');grid(ax)
    c.embed_plot(fig,10,494,342,222,'blocks');plt.close(fig)
    fig,ax=base_plot(325,222,(.21,.23,.72,.67));r=data['repeats'].groupby(['repeat','method']).absolute_error.mean().unstack()
    offsets=np.linspace(-.10,.10,len(r))
    for offset,(_,row) in zip(offsets,r.iterrows()):
        ax.plot([offset,1+offset],[row.A_structure_only,row.F_full_solvai],color=GRID,lw=1.1)
    for x,key,color in [(0,'A_structure_only',MID),(1,'F_full_solvai',TEAL)]:
        ax.scatter(x+offsets,r[key],s=26,color=color,edgecolors='white',lw=.5,zorder=3)
        ax.errorbar(x,r[key].mean(),yerr=r[key].std(ddof=1),fmt='_',ms=10,capsize=3,color=INK,zorder=4)
    ax.axhline(pimd,color=ROSE,lw=.8,ls='--')
    ax.text(.04,pimd-.012,'PIMD8',transform=ax.get_yaxis_transform(),
            ha='left',va='top',color=ROSE,fontsize=8)
    ax.set(xlim=(-.25,1.27),ylim=(.177,.333),xticks=[0,1],xticklabels=['Structure','SolvAI'],ylabel='OOF MAE (kcal mol$^{-1}$)');grid(ax,'y')
    c.embed_plot(fig,389,494,325,222,'repeats');plt.close(fig)
    return c.save(OUT/'F3_matched_evidence',preview/'F3_matched_evidence.png' if preview else None)


def modern(data,preview):
    c=Art(690,'Direct molecular learning and response-augmented pretrained representations')
    c.title(10,25,'a','Test alternative molecular representations')
    c.text(20,52,'Direct graph learning',12,BLUE,600)
    c.molecule(18,66,50,47);c.arrow(78,88,104,88,BLUE)
    c.text(114,83,'D-MPNN / CheMeleon',12,INK,600)
    c.text(114,103,'Fit the graph model end to end',11,MID)
    c.line(357,45,357,120,GRID,.8)
    c.text(380,52,'Frozen encoder + matched endpoint',12,TEAL,600)
    c.text(380,85,'MoLFormer\nor Uni-Mol',11.5,INK,600)
    c.lock(469,70,BLUE)
    c.arrow(482,88,512,88,TEAL)
    c.trees(524,74,87,40)
    c.text(626,86,'With / without',11.5,INK)
    c.text(626,103,'responses',11.5,TEAL,600)
    c.line(10,133,710,133,GRID)
    c.title(10,162,'b','Compare every model on the same held-out molecules')
    order=['structure','solvai','computed_only','dmpnn','chemeleon','molformer_structure','molformer_solvai','unimol_structure','unimol_solvai']
    names=['Structural descriptors','SolvAI · 15 responses','Computation-only · 10','D-MPNN from scratch','CheMeleon fine-tuned','MoLFormer + structure','MoLFormer + responses','Uni-Mol + structure','Uni-Mol + responses']
    for i,name in enumerate(names):
        color=TEAL if i in (1,6,8) else BLUE if i==2 else INK
        c.text(18,221.46+i*25.0278,name,11.7,color,600 if i==1 else 400)
    for k,cohort in enumerate(COHORTS):
        x=245+k*156
        c.text(x+63,191,cohort,12,INK,600,'middle')
        fig,ax=base_plot(156,265,(.09,.12,.86,.85))
        vals=[mae(data,m,cohort) for m in order]
        maxx=.68 if k==0 else 2.7
        for i,(model,v) in enumerate(zip(order,vals)):
            color=TEAL if i in (1,6,8) else BLUE if i==2 else MID
            marker='s' if i==2 else 'o'
            ax.plot([0,v],[i,i],color=GRID,lw=1.1)
            ax.scatter(v,i,s=20,color=color,marker=marker,zorder=3)
            ax.text(v+maxx*.04,i,f'{v:.3f}',va='center',fontsize=7.6,color=color)
        ax.set(ylim=(8.5,-.5),xlim=(0,maxx),yticks=[],xticks=[0,.3,.6] if k==0 else [0,1,2],xlabel='MAE (kcal mol$^{-1}$)')
        ax.spines['left'].set_visible(False);grid(ax)
        c.embed_plot(fig,x,197,156,265,f'scores{k}');plt.close(fig)
    c.title(10,498,'c','Add responses to each fixed representation')
    contrasts=[('solvai','structure'),('molformer_solvai','molformer_structure'),('unimol_solvai','unimol_structure')]
    for i,label in enumerate(['Structure','MoLFormer','Uni-Mol']):c.text(18,561.26+i*34.486,label,11.8,INK)
    for k,cohort in enumerate(COHORTS):
        x=174+k*177;c.text(x+75,529,cohort,12,INK,600,'middle')
        fig,ax=base_plot(172,137,(.10,.25,.84,.73))
        for i,(a,b) in enumerate(contrasts):
            frame=data['pairs'];r=frame.loc[frame.cohort.eq(cohort)&frame.candidate.eq(a)&frame.reference.eq(b)].iloc[0]
            ax.errorbar(r.delta,i,xerr=[[r.delta-r.low],[r.high-r.delta]],fmt='o',ms=4,color=TEAL,capsize=2.5,lw=1.1)
        ax.axvline(0,color=MID,lw=.7);ax.set(ylim=(2.45,-.45),xlim=(-1,.04),yticks=[],xticks=[-1,-.5,0]);ax.spines['left'].set_visible(False);grid(ax)
        c.embed_plot(fig,x,539,172,137,f'augmentation{k}');plt.close(fig)
    c.text(361,672,'ΔMAE after response augmentation (kcal mol⁻¹); negative favors responses',11.5,TEAL,400,'middle')
    return c.save(OUT/'F4_modern_comparisons',preview/'F4_modern_comparisons.png' if preview else None)


def transfer(data,preview):
    c=Art(639,'Solvent-response transfer across chemical separation and external chemistry')
    c.text(10,22,'Structure only',12,MID,600);c.circle(106,18,3.5,MID)
    c.text(133,22,'SolvAI',12,TEAL,600);c.circle(182,18,3.5,TEAL)
    c.title(10,62,'a','Chemical groups')
    c.title(256,62,'b','Nearest neighbors')
    c.title(504,62,'c','No ARROW labels')
    sep=data['separation']
    fig,ax=base_plot(242,224,(.36,.24,.61,.70))
    for i,regime in enumerate(['global_butina_0_70','global_scaffold','global_family']):
        v=sep.loc[sep.regime.eq(regime)].set_index('method').mae
        a,b=v.A_structure_only,v.F_full_solvai
        ax.plot([a,b],[i,i],color=GRID,lw=2)
        ax.scatter(a,i,s=31,color=MID);ax.scatter(b,i,s=31,color=TEAL)
    ax.set(yticks=range(3),yticklabels=['Clusters','Scaffolds','Families'],ylim=(2.5,-.5),xlim=(0,1.4),xticks=[0,.5,1],xlabel='MAE (kcal mol$^{-1}$)');grid(ax)
    c.embed_plot(fig,4,77,242,224,'groups');plt.close(fig)
    fig,ax=base_plot(240,224,(.24,.24,.70,.70))
    for method,color in [('A_structure_only',MID),('F_full_solvai',TEAL)]:
        xs=[.5,.6,.7,.8];ys=[sep.loc[sep.regime.eq(f'global_nn_{x:.2f}')&sep.method.eq(method),'mae'].item() for x in xs]
        ax.plot(xs,ys,color=color,marker='o',ms=4,lw=1.2)
    ax.set(ylim=(.18,.39),xticks=xs,xlabel='Exclusion threshold',ylabel='MAE (kcal mol$^{-1}$)');grid(ax,'y')
    c.embed_plot(fig,255,77,240,224,'neighbors');plt.close(fig)
    fig,ax=base_plot(217,224,(.23,.24,.71,.70))
    z=data['zero'];a=z.loc[z.method.eq('A_structure_only')].sort_values('molecule_id');b=z.loc[z.method.eq('F_full_solvai')].sort_values('molecule_id')
    assert list(a.molecule_id)==list(b.molecule_id)
    delta=np.sort(b.absolute_error.to_numpy()-a.absolute_error.to_numpy())
    ax.axhline(0,color=MID,lw=.7);ax.scatter(range(1,86),delta,color=INK,s=8)
    ax.axhline(delta.mean(),color=INK,lw=.9,ls='--')
    ax.set(xlabel='Ranked molecule',ylabel='Δerror (kcal mol$^{-1}$)',xticks=[1,40,85]);grid(ax,'y')
    c.embed_plot(fig,500,77,217,224,'zeroarrow');plt.close(fig)
    c.text(604,311,f'Mean change: {delta.mean():+.3f}',11.3,INK,600,'middle')
    c.text(604,329,'Δerror: SolvAI − structure',10.5,MID,400,'middle')
    c.line(10,335,710,335,GRID)
    c.title(10,367,'d','External matched comparison')
    c.title(385,367,'e','Relate error to chemical distance')
    fig,ax=base_plot(346,199,(.22,.22,.74,.69))
    for i,cohort in enumerate(['External-220','Strict-97']):
        a,b=mae(data,'structure',cohort),mae(data,'solvai',cohort)
        ax.plot([b,a],[i,i],color=GRID,lw=2)
        ax.scatter(a,i,color=MID,s=38);ax.scatter(b,i,color=TEAL,s=38)
        ax.text(a,i-.18,f'{a:.3f}',ha='center',fontsize=8.1,color=MID)
        ax.text(b,i+.24,f'{b:.3f}',ha='center',fontsize=8.1,color=TEAL)
    ax.set(yticks=[0,1],yticklabels=['220','97'],ylim=(1.5,-.5),xlim=(.9,2.35),xlabel='MAE (kcal mol$^{-1}$)');grid(ax)
    c.embed_plot(fig,10,384,346,199,'external');plt.close(fig)
    for i,key in enumerate(['endpoint_disjoint','strict_response_source_disjoint']):
        r=data['metrics']['external_validation'][key];p=r['paired_difference'];lo,hi=p['ci95']
        c.text(22,598+i*20,f"{r['n']}: ΔMAE {p['mean']:+.3f} [{lo:+.3f}, {hi:+.3f}]",11.3,TEAL)
    fig,ax=base_plot(326,224,(.19,.25,.78,.69));s=data['strata'];xx=np.arange(4)
    ax.bar(xx-.19,s.structure_mae,.35,color=MID);ax.bar(xx+.19,s.solvai_mae,.35,color=TEAL)
    ax.set(xticks=xx,xticklabels=[f'{a:.1f}–{b:.1f}\nn={n}' for a,b,n in zip(s.lower,s.upper,s.n)],ylim=(0,4.6),ylabel='MAE (kcal mol$^{-1}$)');grid(ax,'y')
    c.embed_plot(fig,385,384,326,224,'similarity');plt.close(fig)
    c.text(552,630,'Nearest-training Morgan similarity',11.4,MID,400,'middle')
    return c.save(OUT/'F5_generalization',preview/'F5_generalization.png' if preview else None)


def composition(data,preview):
    c=Art(567,'Response-source composition changes the preferred model across chemistry')
    c.title(10,25,'a','Keep the learner fixed, change the response sources')
    names=['COSMOtherm','Abraham','OpenFF','GBn2','SMD(water)','ConfSolv']
    for i,(name,n) in enumerate(zip(names,COUNTS)):
        x=145+i*93
        c.text(x+37,66,name,11.4,COLORS[i],600,'middle')
        cell=(74-(n-1)*3)/n
        for j in range(n):
            c.rect(x+j*(cell+3),89,cell,14,COLORS[i])
        if i in (2,3):c.text(x+37,100,'physical + δ',10,'white',600,'middle')
        if i==1:
            c.text(x+37,145,'Omitted',11.5,MID,400,'middle')
        else:
            for j in range(n):c.rect(x+j*(cell+3),132,cell,14,COLORS[i])
            if i in (2,3):c.text(x+37,143,'physical',10,'white',600,'middle')
    c.text(10,94,'Full model',12.5,TEAL,600);c.text(10,111,'15 responses',11.5,TEAL)
    c.text(10,138,'Computation-only',12.5,BLUE,600);c.text(10,155,'10 responses',11.5,BLUE)
    c.text(10,188,'Remove five Abraham values and two residual corrections (δ). Retain the two physical energies.',11.3,MID)
    c.text(10,208,'Both models learn hydration from the same experimental labels.',11.5,MID)
    c.line(10,226,710,226,GRID)
    c.title(10,257,'b','Compare absolute errors')
    c.title(389,257,'c','Resolve the change in ordering')
    fig,ax=base_plot(349,261,(.18,.23,.77,.72));xx=np.arange(3)
    for shift,model,color,label,marker in [(-.16,'structure',MID,'Structure','o'),(0,'solvai',TEAL,'Full 15','o'),(.16,'computed_only',BLUE,'Computation-only','s')]:
        vals=[mae(data,model,cohort) for cohort in COHORTS]
        ax.scatter(xx+shift,vals,color=color,s=37,marker=marker,label=label,zorder=3)
    ax.set(xticks=xx,xticklabels=['ARROW\n85','External\n220','Strict\n97'],ylim=(0,2.45),xlim=(-.5,2.5),ylabel='MAE (kcal mol$^{-1}$)');grid(ax,'y')
    ax.legend(loc='upper left',fontsize=8,frameon=False,handlelength=.9,labelspacing=.4)
    c.embed_plot(fig,10,282,349,261,'absolute');plt.close(fig)
    fig,ax=base_plot(326,261,(.31,.23,.65,.72))
    for i,cohort in enumerate(COHORTS):
        f=data['pairs'];r=f.loc[f.cohort.eq(cohort)&f.candidate.eq('solvai')&f.reference.eq('computed_only')].iloc[0]
        color=TEAL if r.delta<0 else BLUE
        ax.errorbar(r.delta,i,xerr=[[r.delta-r.low],[r.high-r.delta]],fmt='o',color=color,ms=5,capsize=3,lw=1.25)
        ax.text(r.delta,i+.25,f'{r.delta:+.3f}',ha='center',color=color,fontsize=8.3)
    ax.axvline(0,color=MID,lw=.8)
    ax.set(yticks=range(3),yticklabels=COHORTS,ylim=(2.6,-.55),xlim=(-.09,.20),xticks=[-.05,0,.05,.1,.15,.2],xlabel='Full 15 − computation-only\nΔMAE (kcal mol$^{-1}$)');grid(ax)
    c.embed_plot(fig,385,282,326,261,'composition');plt.close(fig)
    c.text(431,561,'← Full 15 better',11.3,TEAL)
    c.text(707,561,'Computation-only better →',11.3,BLUE,400,'end')
    return c.save(OUT/'F6_source_composition',preview/'F6_source_composition.png' if preview else None)


def protocol(data,preview):
    c=Art(459,'Inner model selection and fresh graph-model refitting')
    c.title(10,25,'a','Keep test labels outside fitting and selection')
    c.text(16,63,'Outer-training labels',14,INK,600)
    c.rect(16,85,170,21,'#D6E6F2');c.rect(186,85,30,21,'#F0D9B6')
    c.text(101,100,'85% inner fit',11.5,BLUE,600,'middle')
    c.text(201,126,'15%',11.5,AMBER,600,'middle')
    c.text(201,143,'validation',10.7,AMBER,400,'middle')
    c.arrow(230,96,271,96,BLUE)
    c.text(295,84,'Select epoch count',13,BLUE,600)
    c.text(295,107,'Lowest validation MAE',11.5,INK)
    c.text(295,128,'Cap 80 · patience 12',11.5,MID)
    c.arrow(458,96,501,96,BLUE)
    c.text(526,84,'Initialize afresh',13,TEAL,600)
    c.text(526,107,'Refit on all outer-training labels',11.2,INK)
    c.text(526,128,'Use the selected epoch count',11.2,MID)
    c.text(16,178,'The 15% inner split is selected separately within each label source.',11.5,MID)
    c.line(10,201,710,201,GRID)
    c.title(10,230,'b','Average independent seeds before scoring')
    for i,seed in enumerate([11,29,47]):
        x=30+i*124
        c.text(x+41,271,f'Seed {seed}',12,INK,600,'middle')
        c.molecule(x+9,279,69,50,hydrogens=False)
        c.path(f'M{x+43},327 C{x+43},360 424,347 451,347',BLUE,1)
    c.text(248,390,'Repeat selection and refit for each seed',11.4,MID,400,'middle')
    c.text(474,310,'Average',13,TEAL,600)
    c.text(474,330,'predictions',13,TEAL,600)
    c.arrow(452,347,530,347,TEAL)
    c.arrow(544,347,585,347,TEAL)
    c.text(607,317,'Score',14,INK,600)
    c.text(607,339,'Held-out',11.7,INK)
    c.text(607,356,'molecules',11.7,INK)
    c.text(588,423,'Held-out labels',12,AMBER,600)
    c.arrow(631,404,631,367,AMBER)
    c.text(16,441,'Five ARROW folds + one external fit; three seeds × two direct graph models = 36 runs.',11.5,MID)
    return c.save(SI/'Supp_Fig6_evaluation_protocol',preview/'Supp_Fig6_evaluation_protocol.png' if preview else None)


def conformers(data,preview):
    c=Art(553,'Conformer-resolved sources and a molecular equilibrium hydration target')
    c.title(10,25,'a','Many conformers, one equilibrium hydration target')
    for i,m in enumerate(c.conformers):
        c.molecule(20+i*128,51,108,71,mol=m,angle=-15+i*25)
        c.text(74+i*128,145,f'Conformer {i+1}',11.5,MID,400,'middle')
    c.arrow(413,93,464,93,TEAL)
    c.text(581,86,'One molecular property',13,INK,600,'middle')
    c.text(559,122,'ΔG',26,TEAL,600);c.text(603,129,'hyd',11,TEAL)
    c.line(10,167,710,167,GRID)
    c.title(10,195,'b','Summarize the supplied ensembles')
    for x,label,color,levels in [(42,'Gas',MID,[250,296,325]),(237,'Water',BLUE,[297,250,325])]:
        c.text(x+43,227,label,13,color,600,'middle')
        for j,yy in enumerate(levels):
            c.line(x,yy,x+55,yy,color,2)
            c.text(x+64,yy+4,f'g{j+1}',12,color)
    for j,(a,b) in enumerate(zip([250,296,325],[297,250,325])):
        c.path(f'M138,{a} C172,{a} 192,{b} 225,{b}',GRID,.8)
    c.text(180,354,'Schematic levels; not measured source energies',11,MID,400,'middle')
    fig=plt.figure(figsize=(79/25.4,40/25.4));fig.text(.03,.77,r'$g_i=G_i-\min_jG_j$',fontsize=12,color=INK)
    fig.text(.03,.41,r'$C=-RT\ln\sum_i e^{-g_i/RT}$',fontsize=12,color=INK)
    c.embed_plot(fig,380,230,316,160,'equations');plt.close(fig)
    c.title(10,403,'c','Predict six summaries from molecular structure')
    mathlabels=[r'$C_{\mathrm{gas}}$',r'$C_{\mathrm{water}}$',r'$\Delta G_{\mathrm{ensemble}}-s_{i_0}$',r'$\mathrm{SD}(s_i)$',r'$\mathrm{mean}(\Delta g_i)$',r'$\mathrm{SD}(\Delta g_i)$']
    desc=['Gas ensemble','Water ensemble','Hydration correction','Solvation spread','Mean water response','Response spread']
    for i,(mathtext,label) in enumerate(zip(mathlabels,desc)):
        x=10+i*118
        c.rect(x,425,110,46,'#F4F0F7',radius=1)
        fig=plt.figure(figsize=(27.5/25.4,11.5/25.4));fig.text(.5,.42,mathtext,fontsize=8.6,color=PURPLE,ha='center')
        c.embed_plot(fig,x,425,110,46,f'summary{i}');plt.close(fig)
        c.text(x+55,489,label,10.5,MID,400,'middle')
    fig=plt.figure(figsize=(175/25.4,12/25.4))
    fig.text(0,.62,r'$s_i$: conformer solvation energy; $i_0$: lowest-gas-energy conformer.',fontsize=8.1,color=MID)
    fig.text(0,.15,r'$\Delta g_i=g_i^{\rm water}-g_i^{\rm gas}$. Each phase uses its own energy minimum.',fontsize=8.1,color=MID)
    c.embed_plot(fig,10,503,700,48,'definitions');plt.close(fig)
    return c.save(SI/'Supp_Fig7_conformer_targets',preview/'Supp_Fig7_conformer_targets.png' if preview else None)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--preview-dir',type=Path)
    functions=[overview,learning,matched,modern,transfer,composition,protocol,conformers]
    parser.add_argument('--only',nargs='+',choices=[f.__name__ for f in functions])
    args=parser.parse_args()
    data=load_data()
    for function in functions:
        if args.only and function.__name__ not in args.only:continue
        result=function(data,args.preview_dir)
        print(f'Built locally on {platform.system()} {platform.machine()}: {result}',flush=True)
    inputs=['results/paper_metrics.json',
        'results/confirmatory/standardized_exclusion_endpoint_predictions.parquet',
        'results/confirmatory/confirmatory_paired_comparisons.csv',
        'results/confirmatory/standardized_exclusion_global_separation_metrics.csv',
        'results/endpoint_models_20261006/corrected_original_metrics.csv',
        'results/endpoint_models_20261006/corrected_original_comparisons.csv',
        'results/journal_20261002/similarity_strata.csv',
        'figures/source/fig1_assets/nma_openbabel.pdb',
        'figures/source/fig1_assets/dimethoxyethane_selected_conformers.sdf']
    inputs += [str(p.relative_to(ROOT)) for p in sorted((OUT/'teaser_components').glob('*.png'))]
    scripts=['scripts/make_local_journal_figures.py','scripts/journal_vector_art.py','scripts/journal_teaser.py','scripts/journal_style.py']
    manifest=dict(source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        builder_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in scripts},
        build_platform=f'{platform.system()} {platform.machine()}',
        source_descriptor_counts=COUNTS,computed_only_descriptor_count=10,
        embedded_raster_images=True,mixed_figures=['F1_overview','F2_learning_and_evaluation'],molecular_connectivity_checked=True,
        molecular_illustrations='N-methylacetamide; 1,2-dimethoxyethane conformers',
        solvent_shells_and_energy_levels='Schematic, not simulation results',
        main_figures=[p.stem for p in sorted(OUT.glob('F*.pdf'))])
    (OUT/'sources.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
