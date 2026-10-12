#!/usr/bin/env python3
"""Separate finite-upper-limit boundary terms from pinned CRPropa quadrature.

This diagnostic adds no production tail. It uses the same zero-truncated cross
section. Integration by parts has a nonzero endpoint term if the photon
integral above Kmax/(2*gamma) is nonzero, even when sigma(Kmax)=0.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from diagnose_network_photo_rates import PhotonIntegral, interpolate, MPC_CM
from run_network_rate_scan import verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--reference-data',type=Path,required=True)
    p.add_argument('--channels',default='11:4:full,9:4:full,56:26:full')
    p.add_argument('--gammas',default='1e10,1e11,1e12,1258925411794.1663,1318256738556.4102,1e13,1e14')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();verify(a.data,'manifest.json');verify(a.reference_data,'shared_manifest.json')
    sys.path.insert(0,str(a.reference_data/'sources'))
    import interactionRate
    from units import eV
    energy,density=np.loadtxt(a.reference_data/'shared_photons.tsv').T;photons=PhotonIntegral(energy,density)
    class Field:
        def getDensity(self,x,z=0):return interpolate(np.asarray(x)/eV,energy,density)*1e6/eV
    tables={}
    for filename,offset in [('total.tsv',2),('channels.tsv',3)]:
        for line in (a.data/filename).read_text().splitlines():
            if not line or line.startswith('#'):continue
            r=line.split();tables[int(r[0]),int(r[1]),'full' if offset==2 else r[2]]=np.array(list(map(float,r[offset:])))
    gammas=np.array([float(v) for v in a.gammas.split(',')]);rows=[]
    for item in a.channels.split(','):
        mass,charge,code=item.split(':');mass,charge=int(mass),int(charge)
        grid=np.loadtxt(a.data/('grid_external.tsv' if mass<12 else 'grid_talys.tsv'))*1e6;xs=tables[mass,charge,code]
        if xs[-1]!=0:raise ValueError('This check requires the agreed zero-ended grid')
        nodes,weights=leggauss(64);lo,hi=np.log(grid[:-1]),np.log(grid[1:]);x=np.exp(((hi+lo)[:,None]+(hi-lo)[:,None]*nodes)/2).ravel();w=((hi-lo)[:,None]/2*weights).ravel()
        sigma=interpolate(x,grid,xs);F=np.sum(w*x*x*sigma*1e-27)
        independent=np.array([np.sum(w*x*x*sigma*1e-27*photons(x/(2*g)))*MPC_CM/(2*g*g) for g in gammas])
        boundary=F*photons(grid[-1]/(2*gammas))*MPC_CM/(2*gammas*gammas)
        official=interactionRate.calc_rate_eps(grid*eV,xs*1e-31,gammas,Field());curves={}
        for factor in [8,16,32,64]:
            refined=np.exp(np.linspace(np.log(grid[0]),np.log(grid[-1]),(len(grid)-1)*factor+1))
            curves[factor]=interactionRate.calc_rate_eps(refined*eV,interpolate(refined,grid,xs)*1e-31,gammas,Field())
        for i,g in enumerate(gammas):
            rows.append(dict(a=mass,z=charge,channel=code,gamma=float(g),rest_upper_eV=float(grid[-1]),independent=float(independent[i]),
                boundary=float(boundary[i]),boundary_fraction=float(boundary[i]/independent[i]) if independent[i] else None,
                official=float(official[i]),refined_plus_boundary={str(f):float(v[i]+boundary[i]) for f,v in curves.items()},
                corrected_relative={str(f):float((v[i]+boundary[i])/independent[i]-1) if independent[i] else None for f,v in curves.items()}))
    a.output.write_text(json.dumps(dict(status='endpoint/quadrature diagnostic only; production unchanged',rows=rows),indent=2)+'\n')
    print('max corrected factor-64 relative difference:',max(abs(r['corrected_relative']['64']) for r in rows if r['independent']>.1))


if __name__=='__main__':main()
