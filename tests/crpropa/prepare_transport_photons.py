#!/usr/bin/env python3
"""Add native CRPropa photon tables to an immutable shared-reference copy.

Use the pinned official convolution, including its documented endpoint error.
This creates a control, not a corrected production model or acceptance target.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import numpy as np
from diagnose_network_photo_rates import interpolate
from run_network_rate_scan import verify


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emission_ratio(numerator, denominator):
    # Match pinned processEmission: division followed by np.nan_to_num.
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.nan_to_num(numerator / denominator)


def prepare(data, reference, output):
    for name in ['manifest.json','decay_manifest.json','photon_manifest.json']:
        verify(data,name)
    verify(reference,'shared_manifest.json')
    if output.exists(): raise ValueError('Output exists; preserve immutable references')
    sys.path.insert(0,str(reference/'sources'))
    import interactionRate
    from units import eV
    import scipy
    energy,density=np.loadtxt(reference/'shared_photons.tsv').T
    class Field:
        def getDensity(self,x,z=0):
            return interpolate(np.asarray(x)/eV,energy,density)*1e6/eV
    gamma=np.logspace(6,14,201)
    grid=np.loadtxt(data/'grid_talys.tsv')*1e6*eV
    transitions={}
    for line in (data/'pd_photon_transitions.tsv').read_text().splitlines():
        if not line or line.startswith('#'): continue
        r=line.split();key=tuple(map(int,r[:4]))
        transitions[key]=np.array(list(map(float,r[4:])))*1e-31
    reachable={tuple(v) for v in json.loads((data/'manifest.json').read_text())['reachable_states']}
    denominators={};rows=[];saturated=0;zero_denominator=0
    for line in (data/'pd_photons.tsv').read_text().splitlines():
        if not line or line.startswith('#'): continue
        r=line.split();key=tuple(map(int,r[:4]));a,z,ar,zr=key
        if (a,z) not in reachable: continue
        if key not in denominators:
            denominators[key]=interactionRate.calc_rate_eps(grid,transitions[key],gamma,Field())
        numerator=interactionRate.calc_rate_eps(grid,np.array(list(map(float,r[5:])))*1e-31,gamma,Field())
        denominator=denominators[key]
        if np.any(~np.isfinite(numerator)) or np.any(numerator<0) or np.any(~np.isfinite(denominator)) or np.any(denominator<0):
            raise ValueError('Invalid official photon convolution')
        ratio=emission_ratio(numerator,denominator)
        saturated+=int(np.count_nonzero(ratio>1));zero_denominator+=int(np.count_nonzero((numerator>0)&(denominator==0)))
        rows.append(f'{z} {a-z} {zr} {ar-zr} {float(r[4])*1e6:.9g} '+' '.join(f'{v:.9g}' for v in ratio))
    shutil.copytree(reference,output,ignore=shutil.ignore_patterns('__pycache__'))
    target=output/'Photodisintegration/photon_emission_Sha.txt'
    target.write_text('# Pinned official convolution; CRPropa Bernoulli saturation; no energy debit\n'+'\n'.join(rows)+'\n')
    corrected=data/'decay_sources/corrected/data/nuclear_decay.txt'
    if sha(corrected)!=sha(output/'nuclear_decay.txt'): raise ValueError('Corrected decay reference differs')
    d=dict(status='official CRPropa photon control; not physical acceptance',revision='34b10e27375c0e9894d7c984ac15fb96e12975ef',
        reference_manifest_sha256=sha(reference/'shared_manifest.json'),photon_manifest_sha256=sha(data/'photon_manifest.json'),
        script_sha256=sha(Path(__file__)),lines=len(rows),saturated_line_gamma_cells=saturated,positive_over_zero_cells=zero_denominator,
        dependencies=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        payload_sha256={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='shared_manifest.json'})
    # Replace the copied manifest after changing the photon payload.
    (output/'shared_manifest.json').write_text(json.dumps(d,indent=2)+'\n')
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(prepare(a.data.resolve(),a.reference.resolve(),a.output.resolve()),indent=2))
