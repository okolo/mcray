#!/usr/bin/env python3
"""External-server deterministic PD refinement at the initial pilot point."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from run_pd_reference import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--external-server',action='store_true')
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--references',type=Path,required=True,help='Prepared transport-photons directory containing cmb/ and cmb-inoue/')
    a=p.parse_args()
    if not a.external_server or platform.system()!='Linux' or platform.machine()!='x86_64':p.error('Requires an external native Linux amd64 server')
    root=Path(__file__).resolve().parents[2];out=Path(os.environ['MCRAY_RESULTS_DIR'])/'deterministic-pd'
    out.mkdir(exist_ok=False);cases=[]
    for bg in ['cmb','cmb-inoue']:
        for provider,setting in [('I',16),('I',32),('I',64),('M',.01),('M',.005),('M',.0025),('C',0)]:
            cases.append(dict(background=bg,provider=provider,setting=setting,id=f'{bg}-{provider}-{setting}'))
    manifest=dict(status='running',histories=0,cases=cases,driver_sha256=sha(Path(__file__)),completed=[])
    marker=out/'matrix.json'
    def save():marker.write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        subprocess.run([sys.executable,'-m','unittest','discover','-s','tests/deterministic','-v'],cwd=root,check=True)
        for case in cases:
            command=[sys.executable,str(Path(__file__).with_name('run_pd_reference.py')),
                '--data',str(a.data.resolve()),'--reference',str(a.references.resolve()/case['background']),
                '--provider',case['provider'],'--energy-per-nucleon','3.9810717055349857e18',
                '--distance-mpc','1','--photons','On','--path-control','--output',str(out/case['id']),
                '--photon-bin-edges','0,1e14,1e15,1e16,inf']
            if case['provider']=='I':command+=['--order',str(case['setting'])]
            if case['provider']=='M':command+=['--native-probe',str(root/'build/Release/tests/deterministic_pd_export'),'--background-step',str(case['setting'])]
            with (out/(case['id']+'.log')).open('w') as log:subprocess.run(command,cwd=root,check=True,stdout=log,stderr=log)
            manifest['completed'].append(case['id']);save()
        manifest['status']='complete deterministic diagnostics; no acceptance verdict'
    except BaseException as error:
        manifest['status']='failed or incomplete';manifest['error']=str(error);raise
    finally:save()


if __name__=='__main__':main()
