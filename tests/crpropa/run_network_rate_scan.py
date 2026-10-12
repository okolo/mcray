#!/usr/bin/env python3
"""Run a frozen Fe56 network rate matrix; this is not physical acceptance.

Use tests/server/run_tests.py to fetch/build/test/archive on external Linux.
The extended profile is intentionally restricted to external Linux amd64.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

PIN = '34b10e27375c0e9894d7c984ac15fb96e12975ef'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def verify(directory, name):
    manifest = json.loads((directory/name).read_text())
    if manifest['revision'] != PIN: raise ValueError('Wrong pinned data revision')
    for file, expected in manifest['payload_sha256'].items():
        if sha(directory/file) != expected: raise ValueError('Checksum mismatch: '+file)
    return manifest


def read_rates(path, gammas, data):
    """Require every expected full/sum/channel row, retaining failures separately."""
    manifest = json.loads((data/'manifest.json').read_text())
    identities = {tuple(v) for v in manifest['reachable_states'] if v[0]>1}
    channels = {id: {'full', 'sum'} for id in identities}
    for line in (data/'channels.tsv').read_text().splitlines():
        if not line or line.startswith('#'): continue
        row = line.split(); identity = tuple(map(int,row[:2]))
        if identity in channels: channels[identity].add(row[2])
    rows = {}; failures = []
    for line in path.read_text().splitlines():
        if line.startswith('# FAILED'): failures.append(line)
        if not line or line.startswith('#'): continue
        a,z,g,c,r = line.split(); key = (int(a),int(z),float(g),c); rate = float(r)
        if key in rows or not math.isfinite(rate) or rate<0: raise ValueError('Invalid/duplicate rate')
        rows[key] = rate
    expected = {(a,z,g,c) for (a,z),codes in channels.items() for g in gammas for c in codes}
    missing = expected-rows.keys(); unexpected = rows.keys()-expected
    summary = dict(rows=len(rows),isotopes=len(identities),gamma_points=len(gammas),
                   failed_cells=failures,missing_rows=len(missing),unexpected_rows=len(unexpected))
    if failures or missing or unexpected: raise ValueError('Incomplete rate matrix: '+json.dumps(summary))
    return rows, summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--cmb-reference',type=Path,required=True)
    parser.add_argument('--inoue-reference',type=Path,required=True)
    parser.add_argument('--probe',type=Path,default=Path('build/Release/tests/nuclear_decay_probe'))
    parser.add_argument('--output',type=Path,default=None)
    parser.add_argument('--profile',choices=['default','convergence','extended'],default='convergence')
    parser.add_argument('--pilot',action='store_true',help='Nine log-spaced points for runtime/driver checks only')
    parser.add_argument('--external-server',action='store_true')
    args=parser.parse_args()
    if args.external_server and (platform.system()!='Linux' or platform.machine()!='x86_64'):
        parser.error('External execution requires native Linux amd64')
    if args.profile=='extended' and not args.external_server:
        parser.error('Extended profile requires --external-server')
    if args.profile=='extended' and args.pilot:parser.error('Extended profile cannot be reduced to a pilot')
    output = args.output or (Path(os.environ['MCRAY_RESULTS_DIR'])/'network-rates' if 'MCRAY_RESULTS_DIR' in os.environ else None)
    if output is None: parser.error('Supply --output or MCRAY_RESULTS_DIR')
    output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    data=args.data.resolve();probe=args.probe.resolve()
    refs={'cmb':args.cmb_reference.resolve(),'cmb-inoue':args.inoue_reference.resolve()}
    for name in ['manifest.json','decay_manifest.json','photon_manifest.json']: verify(data,name)
    for ref in refs.values(): verify(ref,'shared_manifest.json')
    count=9 if args.pilot else (801 if args.profile=='extended' else 401)
    gammas=[10.**(6+8*i/(count-1)) for i in range(count)]
    (output/'gamma-grid.txt').write_text(''.join(format(g,'.17g')+'\n' for g in gammas))
    settings=[(.01,1e-4)]
    if args.profile!='default':settings.extend([(.005,1e-4),(.005,1e-5)])
    if args.profile=='extended':settings.append((.0025,1e-6))
    orders='8,16,32,64,128' if args.profile=='extended' else '8,16,32'
    manifest=dict(status='running; numerical diagnostic only',started_utc=datetime.now(timezone.utc).isoformat(),
        host=platform.platform(),profile=args.profile,pilot=args.pilot,gamma_grid=gammas,orders=orders,settings=settings,seeds='none; deterministic rates',
        data_revision=PIN,probe_sha256=sha(probe),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        scripts_sha256={p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('diagnose_network_photo_rates.py'),Path(__file__).with_name('summarize_network_rates.py'),Path(__file__).with_name('diagnose_rate_boundary.py')]},
        input_manifests_sha256={str(p):sha(p) for p in [data/'manifest.json',data/'decay_manifest.json',data/'photon_manifest.json',*[r/'shared_manifest.json' for r in refs.values()]]},jobs=[])
    start=time.monotonic()
    def save(): (output/'run-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def run(name,command):
        remaining=None if args.external_server else 1800-(time.monotonic()-start)
        if remaining is not None and remaining<=0:raise RuntimeError('Local 1800-second limit reached; use external server')
        tick=time.monotonic();job=dict(name=name,command=command,status='running');manifest['jobs'].append(job);save()
        try:
            with (output/(name+'.log')).open('w') as log:
                result=subprocess.run(command,stdout=log,stderr=log,timeout=remaining,
                    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
            job['exit_code']=result.returncode;job['status']='passed' if result.returncode==0 else 'failed'
            if result.returncode:raise RuntimeError('Failed job: '+name)
        finally:job['seconds']=time.monotonic()-tick;save()
        print(name,round(job['seconds'],2),flush=True)
    try:
        for test in ['test_network_photo_rates.py','test_network_rate_scan.py','test_summarize_network_rates.py']:
            run(test,[sys.executable,str(Path(__file__).with_name(test))])
        for background,ref in refs.items():
            photon=ref/'shared_photons.tsv'
            for step,eps in settings:
                name=f'{background}-step{step:g}-eps{eps:g}'
                target=output/(name+'.tsv')
                run(name,[str(probe),'--mode','photo-rates','--data',str(data),'--a','56','--charge','26','--network','1',
                    '--photon-table',str(photon),'--gamma-list',','.join(format(g,'.17g') for g in gammas),
                    '--background-step',str(step),'--eps-rel',str(eps),'--output',str(target)])
                _,summary=read_rates(target,gammas,data);manifest['jobs'][-1]['coverage']=summary;save()
            default=output/(f'{background}-step0.01-eps0.0001.tsv')
            run(background+'-references',[sys.executable,str(Path(__file__).with_name('diagnose_network_photo_rates.py')),
                '--cache-photons','--orders',orders,'--data',str(data),'--reference-data',str(ref),
                '--mcray-rates',str(default),'--output',str(output/(background+'-references.json'))])
            command=[sys.executable,str(Path(__file__).with_name('summarize_network_rates.py')),
                '--reference',str(output/(background+'-references.json')),'--default',str(default),
                '--output',str(output/(background+'-summary.json'))]
            if args.profile!='default':
                command.extend(['--refined',str(output/(background+'-step0.005-eps0.0001.tsv')),
                    '--tight',str(output/(background+'-step0.005-eps1e-05.tsv'))])
            if args.profile=='extended':
                command.extend(['--finest',str(output/(background+'-step0.0025-eps1e-06.tsv'))])
            run(background+'-summary',command)
            run(background+'-boundary',[sys.executable,str(Path(__file__).with_name('diagnose_rate_boundary.py')),
                '--data',str(data),'--reference-data',str(ref),'--output',str(output/(background+'-boundary.json'))])
        manifest['status']='complete finite-grid diagnostic; not physical acceptance'
    except BaseException as error:
        manifest['status']='failed or incomplete';manifest['error']=str(error);raise
    finally:
        manifest['seconds']=time.monotonic()-start
        manifest['outputs_sha256']={str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='run-manifest.json'}
        save()


if __name__=='__main__':main()
