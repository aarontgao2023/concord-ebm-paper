"""Rebuild publication figures from complete mechanism A datasets, no new fitting."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ebmv2-mpl')
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'report/v2/figures'
OUT.mkdir(parents=True,exist_ok=True)
with (ROOT/'results/v2_dev_mechanism_a/standard_paired_metrics.csv').open() as f:
    rows=list(csv.DictReader(f))
order=['IID_H0','REF_H0','REF_H0_N4','STAGE_H0']
labels=['IID','Reference','Reference, 4N','Stage stress']
selected=[next(r for r in rows if r['cell']==cell and r['metric']=='distance_to_truth' and r['group_or_pair']=='group_mean') for cell in order]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.labelcolor':'#24313A','text.color':'#24313A'})
fig,axes=plt.subplots(1,2,figsize=(7.6,3.5),gridspec_kw={'width_ratios':[1,1.15]})
x=np.arange(4)
for offset,key,color,label in [(-.16,'original_mean','#C2542D','Original'),(.16,'repaired_mean','#1F6F8B','Repaired')]:
    axes[0].bar(x+offset,[float(r[key]) for r in selected],width=.3,color=color,label=label)
axes[0].set_xticks(x,labels,rotation=22,ha='right')
axes[0].set_ylabel('Mean distance to true ordering')
axes[0].set_ylim(0,.36)
axes[0].legend(frameon=False,ncol=2,loc='upper left')
axes[0].set_title('A  Recovery error',loc='left',fontweight='bold')
y=np.arange(4)[::-1]
d=np.array([float(r['paired_difference']) for r in selected]);lo=np.array([float(r['ci95_low']) for r in selected]);hi=np.array([float(r['ci95_high']) for r in selected])
axes[1].axvline(0,color='#B2BAC1',lw=1)
axes[1].errorbar(d,y,xerr=[d-lo,hi-d],fmt='o',color='#1F6F8B',capsize=3,ms=5)
axes[1].set_yticks(y,labels)
axes[1].set_xlim(-.006,.008)
axes[1].set_xticks([-.004,0,.004,.008],['-0.004','0','0.004','0.008'])
axes[1].set_xlabel('Repaired minus original; paired 95% CI')
axes[1].set_title('B  Change in recovery error',loc='left',fontweight='bold')
axes[1].grid(axis='x',alpha=.15)
fig.suptitle('Mechanism development: 12 independent datasets per cell',y=1.02,fontsize=11)
fig.tight_layout(w_pad=2)
fig.savefig(OUT/'mechanism_recovery_a.pdf',bbox_inches='tight')
fig.savefig(OUT/'mechanism_recovery_a.png',dpi=180,bbox_inches='tight')
plt.close(fig)
