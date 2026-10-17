#!/usr/bin/env python3
"""External PD+decay grid, quadrature, native-integrator and path refinements."""
import argparse,json,os,platform,subprocess,sys
from pathlib import Path
from compare_coupled_refinement import compare
from run_pd_reference import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external-server',action='store_true')
    p.add_argument('--data',type=Path,required=True);p.add_argument('--references',type=Path,required=True)
    p.add_argument('--profile',choices=['initial','full'],default='initial');a=p.parse_args()
    if not a.external_server or platform.system()!='Linux' or platform.machine()!='x86_64':p.error('Requires an external native Linux amd64 server')
    root=Path(__file__).resolve().parents[2];out=Path(os.environ['MCRAY_RESULTS_DIR'])/'deterministic-coupled';out.mkdir(exist_ok=False)
    node_counts=[3,5,9] if a.profile=='initial' else [17,33,65];order=16 if a.profile=='initial' else 32
    cases=[]
    for bg in ['cmb','cmb-inoue']:
        for provider in ['I','M','C']:
            for nodes in node_counts:cases.append(dict(bg=bg,provider=provider,nodes=nodes,order=order,step=.0025,id=f'{bg}-{provider}-n{nodes}'))
        cases.append(dict(bg=bg,provider='I',nodes=node_counts[-1],order=order*2,step=.0025,id=f'{bg}-I-quadrature'))
        cases.append(dict(bg=bg,provider='M',nodes=node_counts[-1],order=order,step=.00125,id=f'{bg}-M-background'))
    manifest=dict(status='running',histories=0,cases=cases,driver_sha256=sha(Path(__file__)),completed=[])
    def save():(out/'matrix.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        subprocess.run([sys.executable,'-m','unittest','discover','-s','tests/deterministic','-v'],cwd=root,check=True)
        for case in cases:
            command=[sys.executable,'tests/deterministic/run_coupled_reference.py','--external-server',
                '--data',str(a.data.resolve()),'--reference',str(a.references.resolve()/case['bg']),
                '--provider',case['provider'],'--energy-per-nucleon','3.9810717055349857e18','--nodes',str(case['nodes']),
                '--order',str(case['order']),'--background-step',str(case['step']),'--path-control','--output',str(out/case['id'])]
            if case['provider']=='M':command+=['--native-probe',str(root/'build/Release/tests/deterministic_pd_export')]
            with (out/(case['id']+'.log')).open('w') as log:subprocess.run(command,cwd=root,check=True,stdout=log,stderr=log)
            manifest['completed'].append(case['id']);save()
        for bg in ['cmb','cmb-inoue']:
            for provider in ['I','M','C']:
                for lo,hi in zip(node_counts,node_counts[1:]):
                    report=compare(out/f'{bg}-{provider}-n{lo}',out/f'{bg}-{provider}-n{hi}')
                    (out/f'{bg}-{provider}-grid-{lo}-{hi}.json').write_text(json.dumps(report,indent=2)+'\n')
            for provider,label in [('I','quadrature'),('M','background')]:
                report=compare(out/f'{bg}-{provider}-n{node_counts[-1]}',out/f'{bg}-{provider}-{label}')
                (out/f'{bg}-{provider}-{label}-refinement.json').write_text(json.dumps(report,indent=2)+'\n')
        manifest['status']='complete PD+decay refinements; no acceptance verdict'
    except BaseException as error:manifest.update(status='failed or incomplete',error=str(error));raise
    finally:save()


if __name__=='__main__':main()
