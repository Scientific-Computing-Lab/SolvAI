#!/usr/bin/env python3
"""Rebuild endpoint-study tables and native vector plots from released outputs."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'results/endpoint_models_20261006'
PAPER = ROOT/'paper'
FAMILIES = ['tree','resnet','dual','ridge_residual','graph','molformer_deterministic',
            'tabm','tabpfn','molformer_fixed_tree']
LABELS = {'tree':'ExtraTrees','resnet':'Residual MLP','dual':'Dual branch',
          'ridge_residual':'Ridge + neural residual','graph':'Message passing',
          'molformer_deterministic':'MoLFormer + neural head','tabm':'TabM',
          'tabpfn':'TabPFN-3.5','molformer_fixed_tree':'MoLFormer + ExtraTrees',
          'capacity_adapted':'Width-adaptive selection','count_adapted':'Count-adaptive selection'}
IDEAS = ['Leaf averages; bounded output',
         'Residual links; external response gain',
         'Separate inputs; improved size transfer',
         'Linear trend + correction; low Strict error',
         'Size-aware graph; external response gain',
         'Frozen sequence encoding; low external error',
         'Shared neural ensemble; low size error',
         'Label-conditioned transformer; low size error',
         'Frozen encoding + trees; input contrast']
COHORTS = ['ARROW-85','External-220','Strict-97']
COLORS = ['#647581','#417FA8','#C78935','#9870AD','#277476','#5173B4','#B05C7C','#008C7A','#8097A5']
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
    'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,
    'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,
    'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
    'text.color':'#18303C','axes.labelcolor':'#18303C','svg.hashsalt':'solvai-endpoints'})

from journal_style import apply as apply_style, MODEL_COLORS, LINE_STYLES, normalize_svg
COLORS = MODEL_COLORS
apply_style(plt)

def read(name):return pd.read_csv(DATA/name, float_precision='round_trip')

def save(fig, name, supplement=False):
    directory=PAPER/('supplementary/figures' if supplement else 'figures/journal')
    for ext in ['svg','pdf']:
        fig.savefig(directory/f'{name}.{ext}',facecolor='white')
    normalize_svg(directory/f'{name}.svg')
    if PREVIEW:
        PREVIEW.mkdir(parents=True,exist_ok=True)
        fig.savefig(PREVIEW/f'{name}.png',dpi=160,facecolor='white')
    plt.close(fig)

def table(path, columns, header, rows):
    lines=[r'\begin{longtable}{@{}'+columns+r'@{}}',r'\toprule',header+r' \\',r'\midrule']
    if path.stem=='endpoint_comparisons':
        lines += [r'\endfirsthead',
                  r'\multicolumn{5}{@{}l}{\emph{Supplementary Table 14 (continued)}} \\',
                  r'\toprule',header+r' \\',r'\midrule']
    lines.append(r'\endhead')
    previous=None
    for row in rows:
        if path.stem=='endpoint_metrics' and row[0]=='Strict-97' and previous!='Strict-97':
            lines.append(r'\pagebreak')
        lines.append(' & '.join(map(str,row))+r' \\')
        previous=row[0]
    lines += [r'\bottomrule\end{longtable}']
    path.write_text('\n'.join(lines)+'\n')

def paired(a,b):
    d=np.asarray(a)-np.asarray(b);rng=np.random.default_rng(20260828)
    bs=np.concatenate([d[rng.integers(0,len(d),(1000,len(d)))].mean(1) for _ in range(100)])
    lo,hi=np.quantile(bs,[.025,.975]);return float(d.mean()),float(lo),float(hi)

def corrected_original():
    old=ROOT/'results/journal_20261002'
    f=pd.read_csv(old/'all_predictions.csv',float_precision='round_trip')
    new=read('combined_predictions.csv');masks=f.model.str.startswith('molformer_')
    f=f.loc[~masks].copy()
    for condition in ['structure','solvai']:
        g=new.loc[new.model.eq('molformer_fixed_tree_'+condition),f.columns].copy()
        g['model']='molformer_'+condition;f=pd.concat([f,g],ignore_index=True)
    assert len(f)==9*305 and not f.duplicated(['model','id']).any()
    f.to_csv(DATA/'corrected_original_predictions.csv',index=False)
    metrics=[];contrasts=[]
    oldpairs=pd.read_csv(old/'paired_comparisons.csv')
    for cohort in COHORTS:
        g=f.loc[f.part.eq('arrow') if cohort=='ARROW-85' else f.part.eq('external')]
        if cohort=='Strict-97':g=g.loc[g.strict]
        errors={}
        for model,h in g.groupby('model'):
            h=h.sort_values('id');e=h.prediction-h.y;errors[model]=np.abs(e.to_numpy())
            metrics.append(dict(cohort=cohort,model=model,n=len(h),mae=np.abs(e).mean(),
                                rmse=np.sqrt(np.mean(e**2)),median_ae=np.median(np.abs(e)),
                                r2=1-float(np.sum(e**2)/np.sum((h.y-h.y.mean())**2))))
        for _,r in oldpairs.loc[oldpairs.cohort.eq(cohort)].iterrows():
            if 'molformer_' in r.candidate or 'molformer_' in r.reference:
                delta,lo,hi=paired(errors[r.candidate],errors[r.reference])
            else:delta,lo,hi=r.delta,r.low,r.high
            contrasts.append(dict(cohort=cohort,candidate=r.candidate,reference=r.reference,
                                  delta=delta,low=lo,high=hi,n=len(g)//9))
    m=pd.DataFrame(metrics);p=pd.DataFrame(contrasts)
    m.to_csv(DATA/'corrected_original_metrics.csv',index=False)
    p.to_csv(DATA/'corrected_original_comparisons.csv',index=False)
    from summarize_journal_baselines import LABELS as names
    rows=[[r.cohort,names[r.model],r.n,f'{r.mae:.4f}',f'{r.rmse:.4f}',f'{r.median_ae:.4f}'] for r in m.itertuples()]
    table(PAPER/'supplementary/tables/journal_baselines.tex','llrrrr',
          r'Set & Model & $N$ & MAE & RMSE & Median AE',rows)
    rows=[[r.cohort,names[r.candidate]+r' $-$ '+names[r.reference],
           f'${r.delta:.4f}$',f'${r.low:.4f}$',f'${r.high:.4f}$'] for r in p.itertuples()]
    table(PAPER/'supplementary/tables/journal_comparisons.tex','llrrr',
          r'Set & Contrast & $\Delta$MAE & 95\% lower & 95\% upper',rows)

def tables_and_plots():
    m=read('combined_metrics.csv');p=read('combined_comparisons.csv')
    size=read('size_diagnostic_metrics.csv');sc=read('size_diagnostic_comparisons.csv')
    def val(family,cohort,condition='solvai',agg='prediction'):
        return float(m.loc[m.model.eq(family+'_'+condition)&m.cohort.eq(cohort)&m.aggregation.eq(agg),'mae'].item())
    def sz(family,condition='solvai'):
        v=size.loc[size.model.eq(family+'_'+condition)&size.source.eq('all'),'mae']
        return float(v.item()) if len(v) else None
    rows=[]
    for family,idea in zip(FAMILIES,IDEAS):
        rows.append([LABELS[family],idea,*[f'{val(family,c):.3f}' for c in COHORTS],
                     f'{sz(family):.3f}' if sz(family) is not None else '--'])
    # Main table is a float-compatible tabular, not a longtable inside a float.
    path=PAPER/'tables/endpoint_summary.tex'
    table(path,r'>{\raggedright\arraybackslash}p{3.85cm}>{\raggedright\arraybackslash}p{5.25cm}rrrr',
          r'Endpoint & Mechanism and observed pattern & ARROW & External & Strict & Size',rows)
    text=path.read_text().replace('longtable','tabular').replace(r'\endhead','')
    path.write_text(text)
    eligible=FAMILIES+['count_adapted','capacity_adapted']
    rows=[]
    for cohort in COHORTS:
        for family in eligible:
            for condition,label in [('structure','No response'),('solvai','Response')]:
                r=m.loc[m.model.eq(family+'_'+condition)&m.cohort.eq(cohort)&m.aggregation.eq('prediction')].iloc[0]
                rows.append([cohort,LABELS[family],label,f'{r.mae:.4f}',f'{r.rmse:.4f}',
                             f'{r.median_ae:.4f}',f'{val(family,cohort,condition,"common3_prediction"):.4f}'])
    table(PAPER/'supplementary/tables/endpoint_metrics.tex','lp{3.8cm}lrrrr',
          r'Set & Endpoint & Input & MAE & RMSE & Median & 3-seed MAE',rows)
    rows=[]
    original=pd.read_csv(ROOT/'results/journal_20261002/paired_comparisons.csv')
    for cohort in COHORTS:
        tree=original.loc[original.cohort.eq(cohort)&original.candidate.eq('solvai')&original.reference.eq('structure')].iloc[0]
        rows.append([cohort,LABELS['tree'],f'${tree.delta:.4f}$',f'${tree.low:.4f}$',f'${tree.high:.4f}$'])
        for family in eligible:
            if family=='tree':continue
            q=p.loc[p.aggregation.eq('prediction')&p.cohort.eq(cohort)&
                    p.candidate.eq(family+'_solvai')&p.reference.eq(family+'_structure')]
            assert len(q)==1,(cohort,family)
            r=q.iloc[0]
            rows.append([cohort,LABELS[family],f'${r.delta:.4f}$',f'${r.low:.4f}$',f'${r.high:.4f}$'])
    table(PAPER/'supplementary/tables/endpoint_comparisons.tex','lp{5.3cm}rrr',
          r'Set & Endpoint & $\Delta$MAE & Lower & Upper',rows)
    rows=[]
    for family in eligible:
        if sz(family) is None:continue
        response=sz(family);control=sz(family,'structure')
        rows.append([LABELS[family],f'{response:.4f}',f'{control:.4f}',f'${response-control:.4f}$'])
    table(PAPER/'supplementary/tables/endpoint_size.tex','p{5.4cm}rrr',
          r'Endpoint & Response MAE & No-response MAE & $\Delta$MAE',rows)
    coverage=read('tabpfn_interval_coverage.csv')
    table(PAPER/'supplementary/tables/endpoint_uncertainty.tex','lrrrrr',
          r'Set & Nominal & Covered & Observed & Mean width & Median width',
          [[r.cohort,f'{r.nominal_coverage*100:.0f}\\%',f'{r.covered}/{r.n}',
            f'{r.observed_coverage*100:.1f}\\%',f'{r.mean_width_kcal_mol:.3f}',
            f'{r.median_width_kcal_mol:.3f}'] for r in coverage.itertuples()])
    # Matched information contrast: no legacy/mixed-cache or selected duplicates.
    fig,axes=plt.subplots(1,3,figsize=(7.09,3.85),sharey=True)
    fig.subplots_adjust(left=.30,right=.985,top=.88,bottom=.18,wspace=.30)
    for j,(ax,cohort) in enumerate(zip(axes,COHORTS)):
        for i,(family,color) in enumerate(zip(FAMILIES,COLORS)):
            q=p.loc[p.candidate.eq(family+'_solvai')&p.reference.eq(family+'_structure')&p.cohort.eq(cohort)&p.aggregation.eq('prediction')]
            if q.empty:
                old=pd.read_csv(ROOT/'results/journal_20261002/paired_comparisons.csv')
                q=old.loc[old.candidate.eq('solvai')&old.reference.eq('structure')&old.cohort.eq(cohort)]
            r=q.iloc[0]
            ax.axhline(i,color='#EDF1F3',lw=.55,zorder=0)
            ax.errorbar(r.delta,i,xerr=[[r.delta-r.low],[r.high-r.delta]],fmt='o',ms=4.6,color=color,lw=1.45,capsize=0)
        ax.axvline(0,color='#6A7984',lw=.8);ax.grid(axis='x',color='#E0E7EA',lw=.5)
        ax.set_title('abc'[j]+'  '+cohort,loc='left',fontsize=10)
        ax.set_xlabel(r'$\Delta$MAE');ax.set_yticks(range(len(FAMILIES)));ax.tick_params(axis='y',length=0)
        ax.set_yticklabels([LABELS[x] for x in FAMILIES],fontsize=8.2)
        ax.set_ylim(len(FAMILIES)-.5,-.5);ax.spines['left'].set_visible(False)
    fig.text(.63,.035,'kcal mol$^{-1}$; negative favors adding responses',ha='center',fontsize=8)
    save(fig,'F3_endpoint_comparison')
    fig=plt.figure(figsize=(7.09,6.7));gs=fig.add_gridspec(2,2,height_ratios=[1.13,1])
    fig.subplots_adjust(left=.12,right=.97,top=.92,bottom=.19,hspace=.54,wspace=.32)
    ax=fig.add_subplot(gs[0,:]);plotfamilies=FAMILIES[:-1]
    for i,family in enumerate(plotfamilies):
        ax.axhline(i,color='#EDF1F3',lw=.55,zorder=0)
        ax.plot([sz(family),sz(family,'structure')],[i,i],color='#CCD8DD',lw=1.3)
        ax.scatter(sz(family),i,s=25,color=COLORS[i],zorder=3)
        ax.scatter(sz(family,'structure'),i,s=25,facecolors='white',edgecolors=COLORS[i],zorder=3)
    ax.set_title('a  Accuracy on the size-held-out molecules',loc='left')
    ax.set_yticks(range(8));ax.set_yticklabels([LABELS[f] for f in plotfamilies],fontsize=8.2)
    ax.set_ylim(7.6,-.6);ax.set_xlim(1.2,2.35);ax.set_xlabel('MAE (kcal mol$^{-1}$)')
    ax.grid(axis='x',color='#E0E7EA',lw=.5);ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
    # Extra left room for long labels in the top panel only.
    pos=ax.get_position();ax.set_position([.31,pos.y0,.66,pos.height])
    fig.text(.97,.976,'Filled: with responses   |   Open: without responses',ha='right',fontsize=8)
    stress=read('combined_stress_predictions.csv')
    handles=[]
    for k,series in enumerate(['glycine','alanine']):
        ax=fig.add_subplot(gs[1,k])
        for i,family in enumerate(plotfamilies):
            g=stress.loc[stress.model.eq(family+'_solvai')&stress.series.eq(series)].sort_values('n')
            h,=ax.plot(g.n,g.prediction,color=COLORS[i],lw=1.9 if family in ['tree','tabpfn'] else 1.15,ls=LINE_STYLES[i],
                      marker='o' if family in ['tree','tabpfn'] else None,ms=2.6,label=LABELS[family])
            if k==0:handles.append(h)
        ax.set_title(('b' if k==0 else 'c')+'  Capped '+series,loc='left')
        ax.set_xlabel('Number of residues, n');ax.set_xticks([0,3,6,9,12]);ax.set_ylim(-85,0)
        ax.grid(axis='y',color='#E0E7EA',lw=.5)
        if k==0:ax.set_ylabel(r'Predicted $\Delta G_{hyd}$'+'\n(kcal mol$^{-1}$)')
        ax.text(.98,.96,'No experimental reference',transform=ax.transAxes,
                ha='right',va='top',fontsize=7.5,color='#5E737F')
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.51,.008),ncol=3,frameon=False,fontsize=8,handlelength=2.7,columnspacing=1.3)
    save(fig,'F8_size_and_stress')
    fig,axes=plt.subplots(1,2,figsize=(7.09,3.3))
    fig.subplots_adjust(left=.09,right=.98,top=.84,bottom=.23,wspace=.39)
    for k,cohort in enumerate(['External-220','Strict-97']):
        g=coverage.loc[coverage.cohort.eq(cohort)]
        axes[0].plot(g.nominal_coverage*100,g.observed_coverage*100,'o-',color=[COLORS[7],COLORS[0]][k],label=cohort,ms=5,lw=1.6)
        axes[1].plot(g.nominal_coverage*100,g.mean_width_kcal_mol,'o-',color=[COLORS[7],COLORS[0]][k],label=cohort,ms=5,lw=1.6)
    axes[0].plot([50,100],[50,100],ls='--',color='#B7C3CB',lw=1)
    axes[0].set(xlim=(75,100),ylim=(50,100),xticks=[80,95],xlabel='Nominal coverage (%)',ylabel='Observed coverage (%)')
    axes[1].set(xticks=[80,95],xlabel='Nominal coverage (%)',ylabel='Mean width (kcal mol$^{-1}$)')
    axes[0].set_title('a  Native intervals under-cover',loc='left',fontsize=9)
    axes[1].set_title('b  Width does not ensure coverage',loc='left',fontsize=9)
    for ax in axes:ax.grid(color='#E0E7EA',lw=.5)
    axes[1].legend(frameon=False,fontsize=8,loc='upper left')
    save(fig,'Supp_Fig11_uncertainty',True)
    duration=read('training_duration_pairs.csv');width=read('capacity_width_pairs.csv')
    duration=duration.loc[duration.pair_used_in_horizon_choice & ~duration.family.isin(['embedding','embedding_refined'])]
    width=width.loc[width.used_for_decision]
    fig,axes=plt.subplots(2,1,figsize=(7.09,5.8))
    fig.subplots_adjust(left=.31,right=.97,top=.92,bottom=.11,hspace=.60)
    for k,(frame,column) in enumerate([(duration,'long_minus_ordinary'),(width,'large_minus_base')]):
        grouped=frame.groupby(['family','partition'])[column].mean().reset_index()
        fam=list(grouped.family.unique())
        for j,family in enumerate(fam):
            axes[k].axhline(j,color='#EDF1F3',lw=.55,zorder=0)
            v=grouped.loc[grouped.family.eq(family),column].to_numpy()
            axes[k].scatter(v,j+np.linspace(-.16,.16,len(v)),s=27,
                color=COLORS[FAMILIES.index(family)],edgecolors='white',lw=.5,zorder=3)
        axes[k].axvline(0,color='#768995',lw=.8)
        axes[k].axvline(-.005,color='#98AEB8',ls=':',lw=1)
        axes[k].set_yticks(range(len(fam)))
        axes[k].set_yticklabels([LABELS.get(f,f) for f in fam])
        axes[k].set_ylim(len(fam)-.6,-.6)
        axes[k].set_xlabel('Inner-validation MAE change (kcal mol$^{-1}$)')
        axes[k].grid(axis='x',color='#E3E9EC',lw=.6)
        axes[k].spines['left'].set_visible(False);axes[k].tick_params(axis='y',length=0)
    axes[0].set_title('a  Longer cosine schedules',loc='left')
    axes[1].set_title('b  Width 256 to 512',loc='left')
    save(fig,'Supp_Fig12_optimization',True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preview-dir',type=Path)
    PREVIEW=parser.parse_args().preview_dir
    corrected_original();tables_and_plots()
