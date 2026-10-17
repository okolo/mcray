#!/usr/bin/env python3
"""Prepare/solve a deterministic PD-only model, freezing inputs and failures."""
import argparse
from datetime import datetime,timezone
import hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
import numpy,scipy
from pd_coefficients import build
from transport import solve,bin_edges


def sha(p):
    digest=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):digest.update(b)
    return digest.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--provider',choices=['I','M','C'],required=True);p.add_argument('--energy-per-nucleon',type=float,required=True)
    p.add_argument('--distance-mpc',type=float,default=1.);p.add_argument('--photons',choices=['On','Off'],default='On')
    p.add_argument('--order',type=int,default=32);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--photon-bin-edges',default='',help='Diagnostic edges in eV: 0,...,inf; not frozen acceptance bins')
    p.add_argument('--native-probe',type=Path);p.add_argument('--background-step',type=float,default=.0025);p.add_argument('--eps-rel',type=float,default=1e-6)
    p.add_argument('--path-control',action='store_true',help='Also solve with independent DOP853 integration')
    a=p.parse_args()
    if not 0<=a.distance_mpc<1000 or not 0<a.energy_per_nucleon or not numpy.isfinite(a.energy_per_nucleon) or a.order<2:p.error('Invalid path/energy/order')
    if a.provider=='M' and not a.native_probe:p.error('M requires --native-probe')
    edges=[None if s=='inf' else float(s) for s in a.photon_bin_edges.split(',')] if a.photon_bin_edges else []
    bin_edges(edges)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    root=Path(__file__).resolve().parents[2];data=a.data.resolve();reference=a.reference.resolve()
    manifest=dict(status='running; deterministic diagnostic only',started_utc=datetime.now(timezone.utc).isoformat(),host=platform.platform(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        source_status=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=root,text=True).strip(),
        script_sha256={str(f.relative_to(root)):sha(f) for f in [*Path(__file__).parent.glob('*.py'),
            root/'tests/deterministic_pd_export.cpp',root/'tests/CMakeLists.txt',
            root/'tests/crpropa/diagnose_network_photo_rates.py',root/'tests/crpropa/run_network_rate_scan.py']},
        dependencies=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__),settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        inputs_sha256={str(f):sha(f) for f in [data/'manifest.json',reference/'shared_manifest.json',*([data/'photon_manifest.json'] if a.photons=='On' else [])]},histories=0,seeds='none; deterministic means')
    def save(): (out/'run-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        native=None
        if a.provider=='M':
            native=out/'native-coefficients.tsv';probe=a.native_probe.resolve();manifest['native_probe_sha256']=sha(probe)
            command=[str(probe),'--data',str(data),'--photon-table',str(reference/'shared_photons.tsv'),'--energy-per-nucleon',str(a.energy_per_nucleon),'--output',str(native),'--photons',a.photons,'--background-step',str(a.background_step),'--eps-rel',str(a.eps_rel)]
            manifest['native_command']=command;save()
            with (out/'native-export.log').open('w') as log:
                subprocess.run(command,check=True,stdout=log,stderr=log,timeout=600,env={**os.environ,'OMP_NUM_THREADS':'1'})
        model=build(data,reference,a.provider,a.energy_per_nucleon,a.photons=='On',a.order,native,edges)
        (out/'coefficients.json').write_text(json.dumps(model,indent=2,allow_nan=False)+'\n')
        result=solve(model,a.distance_mpc)
        (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        if a.path_control:
            other=solve(model,a.distance_mpc,'ivp');(out/'path-control.json').write_text(json.dumps(other,indent=2,allow_nan=False)+'\n')
            manifest['path_control_absolute_difference']={k:abs(result['rewards'][k]-other['rewards'][k]) for k in result['rewards']}
        manifest['status']='complete PD-only deterministic diagnostic; no acceptance verdict'
        print('states',len(model['states']),'survival',result['initial_state_yield'],'photon means',result['rewards'])
    except BaseException as error:
        manifest['status']='failed or incomplete';manifest['error']=str(error);raise
    finally:
        manifest['seconds']=time.monotonic()-start
        manifest['outputs_sha256']={f.name:sha(f) for f in out.iterdir() if f.is_file() and f.name!='run-manifest.json'};save()


if __name__=='__main__':main()
