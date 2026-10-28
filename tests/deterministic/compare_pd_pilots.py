#!/usr/bin/env python3
"""Compare a deterministic PD or PD+decay result to matching pilot batches.

Standardized residuals are exploratory diagnostics, not simultaneous intervals
or acceptance verdicts. No histories are generated and no means are corrected.
"""
import argparse
import json
import math
from pathlib import Path
import statistics
from run_pd_reference import sha


def compare(run,pilots,background,code,level):
    manifest=json.loads((run/'run-manifest.json').read_text())
    for name,digest in manifest['outputs_sha256'].items():
        if sha(run/name)!=digest:raise ValueError('Changed deterministic output: '+name)
    combined=manifest['status'].startswith('complete PD+decay')
    if not combined and not manifest['status'].startswith('complete PD-only'):raise ValueError('Incomplete deterministic run')
    model=json.loads((run/'coefficients.json').read_text()) if not combined else dict(metadata=dict(
        provider=manifest['settings']['provider'],energy_per_nucleon_eV=manifest['settings']['energy_per_nucleon'],
        photons=manifest['settings']['pd_photons']=='On'))
    decay_photons=combined and manifest['settings']['decay_photons']=='On'
    result=json.loads((run/'result.json').read_text())
    plan=json.loads((pilots/'jobs.json').read_text());samples=[];hashes={}
    for job in plan['jobs']:
        case=job['case']
        if (case['background']!=background or case['process']!=('talys-decay' if combined else 'talys') or
            case['pd_photons']!=model['metadata']['photons'] or case['decay_photons']!=decay_photons or
            case['energy_per_nucleon_eV']!=model['metadata']['energy_per_nucleon_eV'] or
            case['length_mpc']!=result['distance_Mpc'] or job['code']!=code or job['level']!=level):continue
        marker=pilots/(job['id']+'.done.json');done=json.loads(marker.read_text())
        if sha(pilots/(job['id']+'.tsv'))!=done['output_sha256']:raise ValueError('Changed pilot payload')
        if job['histories']<=0:raise ValueError('Invalid pilot size')
        samples.append((job,done['observables']));hashes[marker.name]=sha(marker)
    if len(samples)<2:raise ValueError('Need at least two matching complete batches')
    if len({job['histories'] for job,_ in samples})!=1:raise ValueError('Unequal batch sizes')
    expected={'survival':result['initial_state_yield'],
        'yield:photon':result['rewards']['photon_count'],'energy:photon':result['rewards']['photon_energy_eV']}
    if combined:expected.update({key:value for key,value in result['rewards'].items() if key.startswith(('yield:','energy:'))})
    for key,value in result['states'].items():
        a,z=map(int,key[1:].split('Z'));species=key if a>1 else ('proton' if z else 'neutron')
        expected['yield:'+species]=value['mean_count'];expected['energy:'+species]=value['mean_energy_eV']
    rows=[]
    for key,mean in expected.items():
        values=[v.get(key,0.) for _,v in samples];mc=statistics.mean(values);se=statistics.stdev(values)/math.sqrt(len(values))
        rows.append(dict(observable=key,deterministic_mean=mean,pilot_mean=mc,batch_standard_error=se,
            difference=mc-mean,standardized_residual=(mc-mean)/se if se else None,
            status='exploratory batch diagnostic' if se else 'unresolved constant/zero sample'))
    return dict(status=('PD+decay' if combined else 'PD-only')+' pilot diagnostic; no acceptance verdict',provider=model['metadata']['provider'],code=code,
        background=background,level=level,batches=len(samples),histories=sum(j['histories'] for j,_ in samples),
        source_sha256=dict(deterministic_manifest=sha(run/'run-manifest.json'),pilot_jobs=sha(pilots/'jobs.json'),markers=hashes),
        observables=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--pilots',type=Path,required=True)
    p.add_argument('--background',choices=['cmb','cmb-inoue'],required=True)
    p.add_argument('--code',choices=['mcray','crpropa'],required=True);p.add_argument('--level',type=int,default=2)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=compare(a.run,a.pilots,a.background,a.code,a.level)
    with a.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
