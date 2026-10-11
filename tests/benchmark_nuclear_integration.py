#!/usr/bin/env python3
"""Paired single-thread propagation timings for two frozen executable directories.

This is a performance comparison, not statistical physical acceptance. Both
variants must implement the stable knot-split integration and --seed in CRbeam.
Outputs include complete nuclear histories and production CRbeam EM cascades.
No run may exceed the local 1800-second wall limit.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from prepare_nuclear_photon_data import verify


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def observables(path, production):
    result = {}
    files = sorted((path/'z0').iterdir()) if production else [path]
    for file in files:
        if not file.is_file(): continue
        for line in file.read_text().splitlines():
            if not line or line.startswith('#'): continue
            row=line.split()
            if production:
                species=file.name; energy,weight=map(float,row[:2])
                energy*=weight
            else:
                _,_,species,weight,energy=row;weight,energy=float(weight),float(energy)
            if not math.isfinite(energy) or not math.isfinite(weight) or energy<0 or weight<0:
                raise ValueError('Invalid output: '+str(file))
            entry=result.setdefault(species,dict(rows=0,weight=0.,energy_eV=0.))
            entry['rows']+=1;entry['weight']+=weight;entry['energy_eV']+=energy
    if not result: raise ValueError('Empty propagation output')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before-dir',type=Path,required=True)
    parser.add_argument('--after-dir',type=Path,required=True)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--primaries',type=int,default=64)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--seed',type=int,default=2026101100)
    args=parser.parse_args()
    if not 1<=args.primaries<=2000 or not 2<=args.repeats<=10 or not 0<=args.seed<=4294967200:
        parser.error('Invalid bounded benchmark parameters')
    output=args.output.resolve();data=args.data.resolve()
    if output==ROOT or ROOT in output.parents: parser.error('Results must be outside the public repository')
    verify(data)
    variants={label:p.resolve() for label,p in [('before',args.before_dir),('after',args.after_dir)]}
    binaries={label:{name:digest(folder/name) for name in ('CRbeam','nuclear_scenario')} for label,folder in variants.items()}
    output.mkdir(parents=True,exist_ok=False)
    cases=[('nuclear-cmb','cmb',False),('nuclear-cmb-inoue','cmb-inoue',False),
           ('crbeam-cmb','0',True),('crbeam-cmb-inoue','3',True)]
    plan=dict(status='running; performance diagnostic only',host=dict(platform=platform.platform(),machine=platform.machine()),
        binaries_sha256=binaries,script_sha256=digest(Path(__file__)),threads=1,primaries=args.primaries,repeats=args.repeats,
        data_sha256={str(p.relative_to(data)):digest(p) for p in sorted(data.rglob('*')) if p.is_file() and '__pycache__' not in p.parts},jobs=[])
    start=time.monotonic()
    def save(): (output/'benchmark.json').write_text(json.dumps(plan,indent=2)+'\n')
    def run(label,case,background,production,repeat,warmup=False):
        seed=args.seed+max(repeat,0);count=min(8,args.primaries) if warmup else args.primaries
        stem=f'{case}-{label}-'+('warmup' if warmup else f'r{repeat}')
        result_path=output/(stem if production else stem+'.tsv')
        if production:
            command=[str(variants[label]/'CRbeam'),'--particle','nucleus','--A','56','--Z','26',
                '--photodisintegration','TALYS','--nuclear-decay','On','--photodisintegration-photons','On','--decay-photons','On',
                '--talys-data',str(data),'--z','2.5e-5','--emax','5.6e20','--emin','1e11','--nparticles',str(count),
                '--background',background,'--seed',str(seed),'--thinning','0','--output',str(result_path)]
        else:
            command=[str(variants[label]/'nuclear_scenario'),'--process','talys-decay','--talys-data',str(data),
                '--a','56','--charge','26','--background',background,'--photodisintegration-photons','On','--decay-photons','On',
                '--energy-per-nucleon','1e19','--distance-mpc','1','--primaries',str(count),'--batches','1',
                '--seed',str(seed),'--particle-output','--output',str(result_path)]
        remaining=1800-(time.monotonic()-start)
        if remaining<=0: raise RuntimeError('Local benchmark wall limit reached')
        tick=time.perf_counter()
        with (output/(stem+'.log')).open('w') as log:
            result=subprocess.run(command,cwd=ROOT/'bin',env={**os.environ,'OMP_NUM_THREADS':'1'},
                stdout=log,stderr=log,timeout=min(240,remaining))
        elapsed=time.perf_counter()-tick
        job=dict(case=case,variant=label,repeat=repeat,warmup=warmup,seed=seed,command=command,
            seconds=elapsed,exit_code=result.returncode,output=str(result_path))
        plan['jobs'].append(job);save()
        if result.returncode: raise RuntimeError('Failed benchmark job: '+stem)
        job['observables']=observables(result_path,production);save()
        print(stem,round(elapsed,3),flush=True)
    try:
        for case,background,production in cases:
            for label in variants:run(label,case,background,production,-1,True)
            for repeat in range(args.repeats):
                order=('before','after') if repeat%2==0 else ('after','before')
                for label in order:run(label,case,background,production,repeat)
        summaries=[]
        for case,_,_ in cases:
            selected=[j for j in plan['jobs'] if j['case']==case and not j['warmup']]
            medians={label:statistics.median(j['seconds'] for j in selected if j['variant']==label) for label in variants}
            pairs=[]
            for repeat in range(args.repeats):
                pair={j['variant']:j for j in selected if j['repeat']==repeat}
                a,b=(pair[label]['observables'] for label in ('before','after'))
                same_species=a.keys()==b.keys()
                same_counts=same_species and all(a[s]['rows']==b[s]['rows'] and a[s]['weight']==b[s]['weight'] for s in a)
                energy=max((abs(a[s]['energy_eV']/b[s]['energy_eV']-1) for s in a if s in b and b[s]['energy_eV']>0),default=0.)
                pairs.append(dict(repeat=repeat,same_species=same_species,same_counts_weights=same_counts,max_energy_relative_difference=energy))
            summaries.append(dict(case=case,median_seconds=medians,speedup=medians['before']/medians['after'],pairs=pairs))
        plan['summary']=summaries;plan['status']='complete; performance diagnostic, not physical acceptance';save()
        print(json.dumps(summaries,indent=2))
    except BaseException as error:
        plan['status']='failed or incomplete';plan['failure']=str(error);save();raise


if __name__=='__main__':main()
