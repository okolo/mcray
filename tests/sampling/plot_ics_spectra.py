#!/usr/bin/env python3
"""Plot frozen single-interaction ICS bin spectra and simulation/oracle ratios.

Requires optional NumPy/Matplotlib. No new Monte Carlo or rebinning is performed.
The spectrum ordinate is the bin-average probability density per ln(energy).
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullFormatter

ENERGIES=(1e13,1e14,1e15,1e16,1e18,1e20)
STYLES={'old':dict(color='#c45b27',marker='s',label='mcray: historical NR'),
        'current':dict(color='#156ca4',marker='o',label='mcray: current')}


def series(reference,rows,case,view,gap=False):
    result=[]
    for ref in reference[case].values():
        if ref['view']!=view or ref['kind']!='bin':continue
        lo,hi=ref['lo'],ref['hi']
        if gap:lo,hi=1-hi,1-lo
        # Log density is undefined in a bin touching zero. These zero-edge
        # bins have no physical probability for the frozen finite-s setup.
        if lo<=0 or hi<=lo:continue
        mu=ref['mean']
        if mu<=0:continue
        width=np.log(hi/lo)
        values={'id':ref['id'],'lo':lo,'hi':hi,'center':np.sqrt(lo*hi),
                'theory':mu/width,'reference':mu,'target':ref['target'],
                'reference_error':ref['mean_error']/width}
        for mode in STYLES:
            r=rows[case,ref['id'],mode]
            values[mode]={'density':r['mean']/width,
                         'lower':r['mean_ci'][0]/width,'upper':r['mean_ci'][1]/width,
                         'N':r['N'],'count':r['counts'][1],
                         'ratio':r['mean']/mu,
                         'ratio_lower':r['mean_ci'][0]/mu,
                         'ratio_upper':r['mean_ci'][1]/mu}
        result.append(values)
    return sorted(result,key=lambda r:r['lo'])


def draw(ax,ratio,data,energy,title,xlabel,limits=None,fraction=False):
    scale=1 if fraction else energy
    centers=np.array([r['center']*scale for r in data])
    xlo=np.array([r['lo']*scale for r in data]);xhi=np.array([r['hi']*scale for r in data])
    theory=np.array([r['theory'] for r in data])
    # Draw disconnected horizontal bins: never interpolate an apparent smooth
    # spectrum across unequal widths or empty bins.
    ax.hlines(theory,xlo,xhi,color='#252525',lw=2,label='Theory: integrated Eq. (12)')
    ymax=max(theory)*4;ymin=max(max(theory)*1e-9,1e-12)
    for mode,style in STYLES.items():
        yy=np.array([r[mode]['density'] for r in data]);lower=np.array([r[mode]['lower'] for r in data]);upper=np.array([r[mode]['upper'] for r in data])
        positive=yy>0
        ax.errorbar(centers[positive],yy[positive],yerr=np.maximum(0,np.array([yy-lower,upper-yy])[:,positive]),ls='none',ms=3.6,capsize=1.8,**style)
        zero=(~positive)&(upper>=ymin)&(upper<=ymax)
        ax.scatter(centers[zero],upper[zero],marker='v',color=style['color'],facecolors='none',s=27)
        # Ratios are displayed where the expected count is sufficient to be
        # informative; this display restriction never changes inference.
        useful=np.array([r['reference']*r[mode]['N']>=20 for r in data])
        rr=np.array([r[mode]['ratio'] for r in data]);rl=np.array([r[mode]['ratio_lower'] for r in data]);ru=np.array([r[mode]['ratio_upper'] for r in data])
        ratio.errorbar(centers[useful],rr[useful],yerr=np.maximum(0,np.array([rr-rl,ru-rr])[:,useful]),ls='none',ms=3,capsize=1.4,color=style['color'],marker=style['marker'])
    ax.set(xscale='log',yscale='log',ylim=(ymin,ymax),title=title,ylabel=r'$dP/d\ln E$ (bin average)')
    ratio.set(xscale='log',ylim=(-.08,3.1),ylabel='MC / theory',xlabel=xlabel)
    ratio.axhline(1,color='#333333',lw=.8);ratio.axhspan(.99,1.01,color='#555555',alpha=.2)
    ratio.axhspan(.9,1.1,color='#999999',alpha=.14)
    if limits is not None:ax.set_xlim(*limits);ratio.set_xlim(*limits)
    else:ax.set_xlim(xlo.min()/1.1,xhi.max()*1.1);ratio.set_xlim(ax.get_xlim())
    ax.tick_params(labelbottom=False)
    for a in (ax,ratio):
        a.xaxis.set_major_locator(LogLocator(base=10,numticks=8))
        a.xaxis.set_major_formatter(LogFormatterMathtext(base=10,labelOnlyBase=True))
        a.xaxis.set_minor_locator(LogLocator(base=10,subs=tuple(range(2,10)),numticks=80))
        a.xaxis.set_minor_formatter(NullFormatter())
        a.grid(alpha=.16,which='major')


def save(fig,stem):
    for extension in ('png','svg','pdf'):
        fig.savefig(stem.with_suffix('.'+extension),dpi=170)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    report=json.loads(a.report.read_text());ref=json.loads(a.reference.read_text())
    rows={(r['case'],r['id'],r['sampler']):r for r in report['spectra']}
    exported=[]
    for view,stem,label in [('gamma','ics_photons',r'Secondary photon energy $E_\gamma$ [eV]'),('recoil','ics_recoil_electrons',r'Recoil electron energy $E_e$ [eV]')]:
        fig=plt.figure(figsize=(13,14.4),layout='constrained');grid=fig.add_gridspec(3,2)
        for i,energy in enumerate(ENERGIES):
            sub=grid[i//2,i%2].subgridspec(2,1,height_ratios=(3,1),hspace=.03)
            ax=fig.add_subplot(sub[0]);ratio=fig.add_subplot(sub[1],sharex=ax)
            case='ICS-'+format(energy,'.0e');data=series(ref,rows,case,view)
            draw(ax,ratio,data,energy,r'$E_0=10^{'+str(int(np.log10(energy)))+r'}$ eV',label)
            if i==0:ax.legend(fontsize=8,loc='best')
            for r in data:
                for mode in STYLES:
                    exported.append({'case':case,'view':view,'bin':r['id'],'sampler':mode,'energy_lo_eV':r['lo']*energy,'energy_hi_eV':r['hi']*energy,'reference_probability':r['reference'],'theory_dP_dlnE':r['theory'],'MC_dP_dlnE':r[mode]['density'],'MC_CI_lower':r[mode]['lower'],'MC_CI_upper':r[mode]['upper'],'N':r[mode]['N'],'count':r[mode]['count']})
        particle='photons' if view=='gamma' else 'recoil electrons'
        fig.suptitle('Single-interaction ICS '+particle+' on CMB (T = 2.7255 K)',fontsize=15)
        fig.supxlabel('Hard photon cutoff: 10¹² eV. Open ▽: simultaneous upper limit for zero events. Ratios shown for expected counts ≥20.',fontsize=9)
        save(fig,a.output_dir/stem)
    fig=plt.figure(figsize=(13,6.8),layout='constrained');grid=fig.add_gridspec(1,2)
    for i,energy in enumerate([1e18,1e20]):
        sub=grid[0,i].subgridspec(2,1,height_ratios=(3,1),hspace=.03)
        ax=fig.add_subplot(sub[0]);ratio=fig.add_subplot(sub[1],sharex=ax)
        data=series(ref,rows,'ICS-'+format(energy,'.0e'),'gamma',gap=True)
        # Changing photon-bin coordinates to q requires a different logarithmic
        # width. series() applies this Jacobian at the bin-probability level.
        draw(ax,ratio,data,energy,r'$E_0=10^{'+str(int(np.log10(energy)))+r'}$ eV',r'$q=1-E_\gamma/E_0=E_e/E_0$',limits=(1e-8,.03),fraction=True)
        ax.set_ylabel(r'$dP/d\ln q$ (transformed photon bins)')
        ax.axvspan(1e-8,1e-4,color='#f7ddd1',zorder=0)
        ax.legend(fontsize=8)
    fig.suptitle('ICS photon endpoint: resolved on a logarithmic energy-gap axis',fontsize=15)
    fig.supxlabel('Shading marks the q<10⁻⁴ tail. Open ▽ denotes an upper limit for zero observations.',fontsize=9)
    save(fig,a.output_dir/'ics_photon_endpoint')
    with (a.output_dir/'ics_spectra.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(exported[0]),delimiter='\t');w.writeheader();w.writerows(exported)

if __name__=='__main__':main()
