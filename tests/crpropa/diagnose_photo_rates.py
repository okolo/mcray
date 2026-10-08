#!/usr/bin/env python3
"""Separate CRPropa rest-energy quadrature and gamma interpolation errors.

This diagnostic never changes production tables. The independent convolution
integrates the agreed log-log cross section against the identical photon table;
the photon integral uses a separate dense logarithmic grid. SciPy is required.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import simpson

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,required=True)
p.add_argument('--reference-data',type=Path,required=True)
p.add_argument('--mcray-rates',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
sys.path.insert(0,str(a.reference_data/'sources'))
import interactionRate
from units import eV
photons=np.loadtxt(a.reference_data/'shared_photons.tsv')
e,d=photons.T

def interpolate(x,grid,values):
    j=np.clip(np.searchsorted(grid,x)-1,0,len(grid)-2)
    f=np.log(x/grid[j])/np.log(grid[j+1]/grid[j])
    positive=(values[j]>0)&(values[j+1]>0)
    left=np.log(values[j],where=values[j]>0,out=np.zeros_like(f))
    right=np.log(values[j+1],where=values[j+1]>0,out=np.zeros_like(f))
    value=np.where(positive,np.exp(left*(1-f)+right*f),values[j]*(1-f)+values[j+1]*f)
    return np.where((x>=grid[0])&(x<=grid[-1]),value,0.)

class Field:
    def getDensity(self,x,z=0):return interpolate(np.asarray(x)/eV,e,d)*1e6/eV

# Independently integrate photons on a dense log grid, preserving the agreed
# log-log interpolation and the linear-in-logE zero-interval convention.
pe=np.geomspace(e[0],e[-1],65537)
pd=interpolate(pe,e,d)
from scipy.integrate import cumulative_trapezoid
integral=-cumulative_trapezoid((pd/pe)[::-1],np.log(pe)[::-1],initial=0)[::-1]
def photon_integral(x):
    return np.interp(np.log(x),np.log(pe),integral,left=integral[0],right=0.)

grid=np.loadtxt(a.data/'grid_talys.tsv')*1e6
row=next(l.split() for l in (a.data/'total.tsv').read_text().splitlines() if l.startswith('56\t26\t') or l.startswith('56 26 '))
xs=np.array(list(map(float,row[2:])))
mc={float(r[2]):float(r[4]) for l in a.mcray_rates.read_text().splitlines() if l and not l.startswith('#') for r in [l.split()] if r[3]=='full'}
gamma=np.array(sorted(mc))
original=interactionRate.calc_rate_eps(grid*eV,xs*1e-31,gamma,Field())
curves={}
for factor in (2,4,8,16):
    refined=np.exp(np.linspace(np.log(grid[0]),np.log(grid[-1]),(len(grid)-1)*factor+1))
    agreed=interpolate(refined,grid,xs)
    curves[str(factor)]=interactionRate.calc_rate_eps(refined*eV,agreed*1e-31,gamma,Field())
    if factor==16:
        independent=np.array([simpson(refined**2*agreed*1e-27*photon_integral(refined/(2*g)),x=np.log(refined))*3.085677581491367e24/(2*g*g) for g in gamma])
rows=[]
for i,g in enumerate(gamma):
    rows.append({'gamma':float(g),'mcray':mc[g],'crpropa_513':float(original[i]),
      'crpropa_refined_agreed_xs':{f:float(v[i]) for f,v in curves.items()},'independent':float(independent[i])})
result={'status':'numerical diagnostic; no production changes or acceptance verdict',
 'inputs_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.mcray_rates,a.data/'total.tsv',a.data/'grid_talys.tsv',a.reference_data/'shared_photons.tsv',a.reference_data/'sources/interactionRate.py',Path(__file__))},
 'photon_integral_nodes':len(pe),'rows':rows}
a.output.write_text(json.dumps(result,indent=2)+'\n')
for threshold in (.1,1.,10.):
    selected=[r for r in rows if r['independent']>threshold]
    print('rate >',threshold,'/Mpc; points',len(selected))
    for key in ('mcray','crpropa_513'):
        print(key,'max relative difference from independent',max(abs(r[key]/r['independent']-1) for r in selected))
    for f in curves:
        print('refinement',f,max(abs(r['crpropa_refined_agreed_xs'][f]/r['independent']-1) for r in selected))
