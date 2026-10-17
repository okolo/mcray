#!/usr/bin/env python3
"""Retain observable-by-observable grid/path refinements, without acceptance."""
import argparse,json,math
from pathlib import Path
from run_pd_reference import sha


def compare(coarse,fine):
    runs=[]
    for directory in [coarse,fine]:
        m=json.loads((directory/'run-manifest.json').read_text())
        if not m['status'].startswith('complete PD+decay'):raise ValueError('Incomplete run')
        for name,digest in m['outputs_sha256'].items():
            if sha(directory/name)!=digest:raise ValueError('Changed output: '+name)
        runs.append((m,json.loads((directory/'result.json').read_text())))
    (a,x),(b,y)=runs
    for key in ['provider','energy_per_nucleon','distance_mpc','pd_photons','decay_photons','photon_bin_edges']:
        if a['settings'][key]!=b['settings'][key]:raise ValueError('Mismatched comparison settings: '+key)
    if a['input_sha256']!=b['input_sha256']:raise ValueError('Changed frozen inputs')
    values=[]
    def add(name,u,v):
        floor=max(a['settings']['atol'],b['settings']['atol'])
        if 'energy' in name:floor*=56*b['settings']['energy_per_nucleon']
        resolved=abs(v)>floor
        values.append(dict(observable=name,coarse=u,fine=v,difference=v-u,
            relative_to_fine=(v-u)/v if resolved else None,
            status='observed refinement; not a certified error bound' if resolved else 'at/below solver coordinate tolerance; no relative precision claim'))
    add('survival',x['initial_state_yield'],y['initial_state_yield'])
    for key in x['rewards']:add(key,x['rewards'][key],y['rewards'][key])
    if set(x['states'])!=set(y['states']):raise ValueError('Changed isotope coverage')
    for key in x['states']:
        for label in ['mean_count','mean_energy_eV']:add(key+':'+label,x['states'][key][label],y['states'][key][label])
    return dict(status='PD+decay numerical refinement diagnostic; no acceptance verdict',
        manifest_sha256=[sha(d/'run-manifest.json') for d in [coarse,fine]],settings=[a['settings'],b['settings']],observables=values)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coarse',type=Path,required=True);p.add_argument('--fine',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=compare(a.coarse,a.fine)
    with a.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
