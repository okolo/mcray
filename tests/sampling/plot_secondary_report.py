#!/usr/bin/env python3
"""Plot selected views of the frozen PP/ICS report (optional Matplotlib)."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True,help='Output stem for PNG and SVG')
    args=parser.parse_args()
    report=json.loads(args.report.read_text());reference=json.loads(args.reference.read_text())
    rows={(r['case'],r['id'],r['sampler']):r for r in report['spectra']}
    styles={'old':dict(color='#bc5426',marker='s',label='Historical NR'),
            'current':dict(color='#146c9d',marker='o',label='Current sampler')}
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.5),layout='constrained')
    def bias(ax,case,view,title):
        bins=sorted((r for r in reference[case].values() if r['view']==view and r['kind']=='bin' and r['mean']>1e-4 and r['lo']>0),key=lambda r:r['lo'])
        x=np.array([np.sqrt(r['lo']*r['hi']) for r in bins])
        for mode,style in styles.items():
            rr=[rows[case,r['id'],mode] for r in bins]
            y=np.array([100*(r['mean']/r['reference']-1) for r in rr])
            lower=np.array([100*r['relative_effect_ci'][0] for r in rr]);upper=np.array([100*r['relative_effect_ci'][1] for r in rr])
            ax.errorbar(x,y,yerr=np.maximum(0,np.array([y-lower,upper-y])),lw=1.3,ms=4,capsize=2,**style)
        ax.axhspan(-1,1,color='#dddddd',label='Core target ±1%',zorder=0)
        ax.axhline(0,color='#777777',lw=.6);ax.set_xscale('log');ax.set_title(title)
        ax.set_xlabel('Secondary energy fraction x');ax.set_ylabel('Bias relative to oracle (%)')
        ax.grid(alpha=.18);ax.legend(fontsize=8)
    bias(axes[0,0],'ICS-1e+15','gamma','ICS, E₀ = 10¹⁵ eV: photon bins')
    bias(axes[0,1],'PP-1e+18','low','Pair production, E₀ = 10¹⁸ eV: lower lepton')
    ax=axes[1,0];case='ICS-1e+20'
    bins=sorted((r for r in reference[case].values() if r['view']=='recoil' and r['kind']=='bin' and r['lo']>0 and r['mean']>1e-12),key=lambda r:r['lo'])
    x=np.array([np.sqrt(r['lo']*r['hi']) for r in bins]);width=np.array([np.log(r['hi']/r['lo']) for r in bins]);p=np.array([r['mean'] for r in bins])
    ax.plot(x,p/width,color='#333333',lw=1.7,label='Numerical oracle')
    for mode,style in styles.items():
        rr=[rows[case,r['id'],mode] for r in bins];y=np.array([r['mean'] for r in rr]);upper=np.array([r['mean_ci'][1] for r in rr]);positive=y>0
        ax.plot(x[positive],(y/width)[positive],ls='--',lw=1.1,ms=4,**style)
        ax.scatter(x[~positive],(upper/width)[~positive],marker='v',color=style['color'],s=32,facecolors='none')
    ax.axvspan(1e-8,1e-4,color='#f7e0d7',zorder=0)
    ax.set(xscale='log',yscale='log',xlabel='Recoil electron fraction Eₑ/E₀',ylabel='Probability density per ln(x)',title='ICS, E₀ = 10²⁰ eV: missing recoil spectrum')
    ax.grid(alpha=.18);ax.legend(fontsize=8);ax.text(.03,.05,'▽: simultaneous upper limit\nfor zero observed events',transform=ax.transAxes,fontsize=8)
    ax=axes[1,1];case='PP-1e+14'
    ids=['low.endpoint.0.000001','low.endpoint.0.000100','low.endpoint.0.010000'];x=np.array([.5-reference[case][key]['lo'] for key in ids])
    ax.plot(x,[reference[case][key]['mean'] for key in ids],color='#333333',lw=1.7,label='Numerical oracle')
    for mode,style in styles.items():
        rr=[rows[case,key,mode] for key in ids];y=np.array([r['mean'] for r in rr]);upper=np.array([r['mean_ci'][1] for r in rr]);positive=y>0
        ax.plot(x[positive],y[positive],ls='--',ms=5,lw=1.1,**style)
        ax.scatter(x[~positive],upper[~positive],marker='v',color=style['color'],s=40,facecolors='none')
    ax.set(xscale='log',yscale='log',xlabel='Distance from equal splitting Δx',ylabel='P(0.5 − x < Δx)',title='Pair production, E₀ = 10¹⁴ eV: equal-split endpoint')
    ax.grid(alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Server spectrum comparison: improved s sampling, shared endpoint defects',fontsize=14)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output.with_suffix('.png'),dpi=180)
    fig.savefig(args.output.with_suffix('.svg'))
    plt.close(fig)

if __name__=='__main__':main()
