#!/usr/bin/env python3
"""Freeze and run pre-acceptance transport/photon pilots on external Linux.

Final families, bins, normalizations, seeds and numerical envelopes must be
approved after these independent pilots. This tool cannot issue acceptance.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

PIN='73fa865b3304113da2059755a83f5939c1039d68'
ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def scenarios(matrix):
    cases=[]
    energies=[10**18.6] if matrix=='initial' else [1e18,10**18.6,1e19]
    distances=[1.] if matrix=='initial' else [1.,10.]
    for bg in ['cmb','cmb-inoue']:
        for energy in energies:
            for length in distances:
                for process,modes in [('talys',[(False,False),(True,False)]),('talys-decay',[(False,False),(True,False),(False,True),(True,True)])]:
                    for pd,decay in modes:
                        cases.append(dict(id=f'{bg}-E{energy:.9g}-L{length:g}-{process}-pd{int(pd)}-dec{int(decay)}',background=bg,
                            energy_per_nucleon_eV=energy,length_mpc=length,process=process,pd_photons=pd,decay_photons=decay))
    return cases


def make_jobs(matrix,profile):
    batches,size=(2,4) if profile=='smoke' else (20,1000)
    result=[]
    for case_index,case in enumerate(scenarios(matrix)):
        for code in ['mcray','crpropa']:
            for level in range(1 if profile=='smoke' else 3):
                for batch in range(batches):
                    # Paired numerical-refinement batches, independent between programs/cases.
                    seed=2026101600+1000000*(code=='crpropa')+100*case_index+batch
                    result.append(dict(case=case,code=code,level=level,batch=batch,histories=size,seed=seed,
                        id=f'{case["id"]}-{code}-n{level}-b{batch:02}'))
    return result


def summarize(path,job):
    histories={};completed=False
    for line in path.read_text().splitlines():
        if line.startswith('# batch=0 primaries='):
            if int(line.split('primaries=')[1])!=job['histories']:raise ValueError('Wrong completion count')
            completed=True
        if not line or line.startswith('#'):continue
        b,h,species,w,e=line.split();b,h=int(b),int(h);w,e=float(w),float(e)
        if b!=0 or not 0<=h<job['histories'] or w!=1 or not math.isfinite(e) or e<0:raise ValueError('Invalid primary/product row')
        row=histories.setdefault(h,Counter())
        if species.startswith('A'):
            a,z=map(int,species[1:].split('Z'))
            if a<2 or not 0<=z<=a:raise ValueError('Invalid nuclear identity')
            row['bound_count']+=w;row['bound_A_sum']+=a*w;row['bound_Z_sum']+=z*w
            row['survival']+=species=='A56Z26'
        else:
            a,z={'proton':(1,1),'neutron':(1,0),'electron':(0,-1),'positron':(0,1),'neutrinoE':(0,0),'neutrinoAE':(0,0),'photon':(0,0)}[species]
        row['baryons']+=a*w;row['charge']+=z*w
        row['yield:'+species]+=w;row['energy:'+species]+=e
        row['energy:all']+=e
        if species=='photon':
            row['photon_histories']=1
            # Exploratory pilot bins only. Preserve individual energies for final bin selection.
            x=e/job['case']['energy_per_nucleon_eV']
            label='underflow' if x<1e-8 else 'overflow' if x>=1 else str(math.floor(math.log10(x)*4))
            row['pilot_photon_bin:'+label]+=w
    if not completed or set(histories)!=set(range(job['histories'])):raise ValueError('Incomplete primary ensemble')
    for row in histories.values():
        if row['baryons']!=56 or row['charge']!=26:raise ValueError('Baryon/charge violation')
    if not job['case']['pd_photons'] and not job['case']['decay_photons']:
        if any(row['yield:photon'] for row in histories.values()):raise ValueError('Photon Off control emitted photons')
    keys=set().union(*(r.keys() for r in histories.values()))
    # Every statistic aggregates descendants inside a complete primary history.
    return {k:sum(r[k] for r in histories.values())/job['histories'] for k in keys}


def reports(output,plan):
    grouped={}
    for job in plan['jobs']:
        marker=output/(job['id']+'.done.json')
        if not marker.exists():continue
        d=json.loads(marker.read_text())
        if sha(output/(job['id']+'.tsv'))!=d['output_sha256']:raise ValueError('Changed completed output')
        grouped.setdefault((job['case']['id'],job['code'],job['level']),{})[job['batch']]=d['observables']
    rows=[];refinements=[];comparisons=[]
    for (case,code,level),batches in sorted(grouped.items()):
        keys=set().union(*(r.keys() for r in batches.values()))
        for key in sorted(keys):
            values=[batches[b].get(key,0.) for b in sorted(batches)]
            se=statistics.stdev(values)/math.sqrt(len(values)) if len(values)>1 else None
            rows.append(dict(case=case,code=code,level=level,observable=key,batches=len(values),mean=statistics.mean(values),batch_standard_error=se,
                uncertainty_status='unresolved zero/constant sample; no zero-width acceptance interval' if se==0 else 'pilot batch estimate'))
        if level>0 and (case,code,level-1) in grouped:
            previous=grouped[case,code,level-1];common=sorted(batches.keys()&previous.keys())
            for key in sorted(keys|set().union(*(r.keys() for r in previous.values()))):
                diffs=[batches[b].get(key,0.)-previous[b].get(key,0.) for b in common]
                refinements.append(dict(case=case,code=code,levels=[level-1,level],observable=key,batches=len(diffs),
                    paired_difference=statistics.mean(diffs),paired_standard_error=statistics.stdev(diffs)/math.sqrt(len(diffs)) if len(diffs)>1 else None))
    for case,code,level in sorted(grouped):
        if code!='mcray' or (case,'crpropa',level) not in grouped:continue
        samples={c:grouped[case,c,level] for c in ['mcray','crpropa']}
        keys=set().union(*(r.keys() for b in samples.values() for r in b.values()))
        for key in sorted(keys):
            values={c:[r.get(key,0.) for r in b.values()] for c,b in samples.items()}
            means={c:statistics.mean(v) for c,v in values.items()}
            se=math.sqrt(sum(statistics.variance(v)/len(v) for v in values.values())) if all(len(v)>1 for v in values.values()) else None
            comparisons.append(dict(case=case,level=level,observable=key,means=means,difference=means['mcray']-means['crpropa'],
                combined_standard_error=se,status='pilot only; zero/constant and rare outputs require separate uncertainty treatment'))
    result=dict(status='pre-acceptance pilots only; no agreement verdict',completed_jobs=sum(len(v) for v in grouped.values()),expected_jobs=len(plan['jobs']),
        observables=rows,numerical_refinements=refinements,program_comparisons=comparisons,
        pending=['approved final family/bins/normalizations','fresh acceptance seeds and fixed sizes','near-zero budgets or inconclusive classification','numerical envelopes with MC uncertainty','simultaneous 99% intervals; 1% integrated/5% spectral budgets'])
    (output/'pilot-report.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--output',type=Path)
    p.add_argument('--matrix',choices=['initial','full'],default='initial');p.add_argument('--profile',choices=['smoke','pilot'],default='pilot')
    p.add_argument('--plan-only',action='store_true');p.add_argument('--external-server',action='store_true');p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if not a.plan_only and (not a.external_server or platform.system()!='Linux' or platform.machine()!='x86_64'):
        p.error('Execution requires an external native Linux amd64 server; use --plan-only locally')
    bundle=a.bundle.resolve()
    inputs=json.loads((bundle/'inputs-manifest.json').read_text())
    for name,h in inputs['payload_sha256'].items():
        if sha(bundle/name)!=h:raise ValueError('Transferred checksum mismatch: '+name)
    output=a.output or (Path(os.environ['MCRAY_RESULTS_DIR'])/'transport-photons' if 'MCRAY_RESULTS_DIR' in os.environ else None)
    if output is None:p.error('Supply --output or MCRAY_RESULTS_DIR')
    output=output.resolve();identity=dict(matrix=a.matrix,profile=a.profile,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scripts_sha256={str(f.relative_to(ROOT)):sha(f) for f in [Path(__file__),ROOT/'tests/nuclear_scenario.cpp',*[ROOT/'tests/crpropa'/name for name in ['decay_probe.cpp','CMakeLists.txt','prepare_transport_photons.py','diagnose_network_photo_rates.py','run_network_rate_scan.py','test_transport_photon_pilots.py','test_transport_photon_tables.py']]]},
        input_manifest_sha256=sha(bundle/'inputs-manifest.json'))
    identity['python']=sys.version
    identity['dependencies']={name:version(name) for name in ['numpy','scipy','GitPython']}
    if a.resume:
        plan=json.loads((output/'jobs.json').read_text())
        if plan['identity']!=identity:raise ValueError('Resume source/tools/inputs/profile changed')
    else:
        output.mkdir(parents=True,exist_ok=False)
        plan=dict(status='pre-acceptance pilots; cannot issue acceptance',identity=identity,crpropa_commit=PIN,jobs=make_jobs(a.matrix,a.profile),started_utc=datetime.now(timezone.utc).isoformat(),host=platform.platform(),cpu_count=os.cpu_count())
        (output/'jobs.json').write_text(json.dumps(plan,indent=2)+'\n')
    if a.plan_only:
        print('Frozen',len(plan['jobs']),'jobs; no calculation started');return
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    def run(label,command,extra=None):
        with (output/(label+'.log')).open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,env={**env,**(extra or {})})
    try:
        for test in ['test_transport_photon_pilots.py','test_transport_photon_tables.py']:
            run(test,[sys.executable,str(ROOT/'tests/crpropa'/test)])
        for bg in ['cmb','cmb-inoue']:
            dest=output/bg
            marker=output/(bg+'-prepared.json')
            if marker.exists():
                from run_network_rate_scan import verify
                verify(dest,'shared_manifest.json')
                if sha(dest/'shared_manifest.json')!=json.loads(marker.read_text())['manifest_sha256']:
                    raise ValueError('Prepared reference manifest changed')
            else:
                if dest.exists():raise ValueError('Incomplete reference preparation; preserve it and start a fresh run')
                run(bg+'-prepare',[sys.executable,str(ROOT/'tests/crpropa/prepare_transport_photons.py'),'--data',str(bundle/'nuclei'),'--reference',str(bundle/bg),'--output',str(dest)])
                marker.write_text(json.dumps(dict(manifest_sha256=sha(dest/'shared_manifest.json')))+'\n')
        reference=output/'reference-source';build=output/'reference-build'
        if not (output/'reference-build.json').exists():
            if reference.exists():raise ValueError('Incomplete reference build; preserve it and start a fresh run')
            run('reference-clone',['git','clone','--no-checkout','https://github.com/CRPropa/CRPropa3.git',str(reference)])
            run('reference-fetch',['git','-C',str(reference),'fetch','origin',PIN]);run('reference-checkout',['git','-C',str(reference),'checkout','--detach','FETCH_HEAD'])
            assert subprocess.check_output(['git','-C',str(reference),'rev-parse','HEAD'],text=True).strip()==PIN
            run('reference-configure',['cmake','-S',str(ROOT/'tests/crpropa'),'-B',str(build),'-DCRPROPA_SOURCE_DIR='+str(reference),'-DBUILD_DOCKER_MCRAY_PROBE=OFF','-DCMAKE_BUILD_TYPE=Release'])
            run('reference-build',['cmake','--build',str(build),'--target','crpropa_decay_probe','--parallel','2'])
            (output/'reference-build.json').write_text(json.dumps(dict(commit=PIN,probe_sha256=sha(build/'crpropa_decay_probe')))+'\n')
        if sha(build/'crpropa_decay_probe')!=json.loads((output/'reference-build.json').read_text())['probe_sha256']:raise ValueError('Reference binary changed')
        native=ROOT/'build/Release/tests/nuclear_scenario'
        native_identity=output/'native-probe.json'
        if native_identity.exists() and json.loads(native_identity.read_text())['sha256']!=sha(native):raise ValueError('Native binary changed')
        native_identity.write_text(json.dumps(dict(sha256=sha(native)))+'\n')
        masses={tuple(map(int,r[:2])):float(r[2]) for line in (bundle/'nuclei/masses.tsv').read_text().splitlines() if line and not line.startswith('#') for r in [line.split()]}
        for i,job in enumerate(plan['jobs']):
            marker=output/(job['id']+'.done.json');target=output/(job['id']+'.tsv')
            if marker.exists():
                if sha(target)!=json.loads(marker.read_text())['output_sha256']:raise ValueError('Completed output changed')
                continue
            c=job['case'];flags=['--photodisintegration-photons','On' if c['pd_photons'] else 'Off','--decay-photons','On' if c['decay_photons'] else 'Off']
            common=['--process',c['process'],'--a','56','--charge','26','--distance-mpc',str(c['length_mpc']),'--primaries',str(job['histories']),'--seed',str(job['seed']),'--output',str(target)]+flags
            if job['code']=='mcray':
                step,eps=[(.01,1e-4),(.005,1e-5),(.0025,1e-6)][job['level']]
                command=[str(native),'--talys-data',str(bundle/'nuclei'),'--photon-table',str(bundle/c['background']/'shared_photons.tsv'),'--particle-output','--pd-photon-yields','CRPropa','--energy-per-nucleon',str(c['energy_per_nucleon_eV']),'--background-step-decades',str(step),'--accuracy',str(eps)]+common
                extra={}
            else:
                command=[str(build/'crpropa_decay_probe'),'--mode','transport','--gamma',str(56*c['energy_per_nucleon_eV']/(masses[56,26]*1e6)),'--step-fraction',str(.03/(3**job['level']))]+common
                extra={'CRPROPA_DATA_PATH':str(output/c['background'])}
            print(f'[{i+1}/{len(plan["jobs"])}] {job["id"]}',flush=True);start=time.monotonic()
            run(job['id'],command,extra)
            values=summarize(target,job)
            marker.write_text(json.dumps(dict(command=command,observables=values,seconds=time.monotonic()-start,output_sha256=sha(target)),indent=2)+'\n')
        plan['status']='complete pre-acceptance pilots; final protocol pending'
    except BaseException as error:
        plan['status']='failed or incomplete';plan['error']=str(error);raise
    finally:
        (output/'jobs.json').write_text(json.dumps(plan,indent=2)+'\n');reports(output,plan)


if __name__=='__main__':main()
