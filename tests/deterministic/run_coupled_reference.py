#!/usr/bin/env python3
"""Freeze and solve PD+decay first moments with an explicit energy grid."""
import argparse
from collections import defaultdict
from datetime import datetime,timezone
import json,os,platform,subprocess,sys,time
from pathlib import Path
import numpy,scipy
from scipy.sparse import save_npz
from pd_coefficients import RawBundle,build,rows,verify
from decay_kernels import DecayBundle,energy_bounds
from coupled_transport import assemble,grids,solve
from run_pd_reference import sha
from transport import bin_edges


def native_decay(path,decay_photons):
    text=path.read_text().splitlines();header=next(s for s in text if s.startswith('# E_per_A_eV='))
    settings=dict(v.split('=') for v in header[2:].split())
    if settings.get('combined')!='On' or settings.get('decay_photons')!=('On' if decay_photons else 'Off'):raise ValueError('Native decay mode mismatch')
    totals={};channels=defaultdict(dict);lines=defaultdict(list)
    for r in rows(path):
        p=tuple(map(int,r[1:3]))
        if r[0]=='T':
            if p in totals:raise ValueError('Duplicate native decay total')
            totals[p]=float(r[3])
        if r[0]=='D':
            if r[3] in channels[p]:raise ValueError('Duplicate native decay channel')
            channels[p][r[3]]=float(r[4])
        if r[0]=='L':lines[(*p,r[3])].append((float(r[4]),float(r[5])))
    return dict(totals=totals,channels=channels,lines=lines)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--provider',choices=['I','M','C'],required=True);p.add_argument('--energy-per-nucleon',type=float,required=True)
    p.add_argument('--distance-mpc',type=float,default=1.);p.add_argument('--pd-photons',choices=['On','Off'],default='On')
    p.add_argument('--decay-photons',choices=['On','Off'],default='On');p.add_argument('--nodes',type=int,default=5)
    p.add_argument('--order',type=int,default=32);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--native-probe',type=Path);p.add_argument('--background-step',type=float,default=.0025);p.add_argument('--eps-rel',type=float,default=1e-6)
    p.add_argument('--photon-bin-edges',default='0,1e14,1e15,1e16,inf')
    p.add_argument('--method',choices=['BDF','Radau'],default='BDF');p.add_argument('--rtol',type=float,default=1e-10);p.add_argument('--atol',type=float,default=1e-16)
    p.add_argument('--path-control',action='store_true');p.add_argument('--external-server',action='store_true')
    a=p.parse_args()
    if not numpy.isfinite(a.energy_per_nucleon) or a.energy_per_nucleon<=0 or not 0<=a.distance_mpc<=10 or a.nodes<2 or a.order<2:p.error('Invalid energy/path/nodes/order')
    if a.provider=='M' and not a.native_probe:p.error('M requires --native-probe')
    if a.external_server and (platform.system()!='Linux' or platform.machine()!='x86_64'):p.error('External mode requires native Linux amd64')
    if not a.external_server and (a.nodes>9 or a.order>32):p.error('Larger refinement grids require --external-server on Linux amd64')
    edges=[None if v=='inf' else float(v) for v in a.photon_bin_edges.split(',')];bin_edges(edges)
    data=a.data.resolve();reference=a.reference.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[2];start=time.monotonic()
    source=[*Path(__file__).parent.glob('*.py'),root/'tests/deterministic_pd_export.cpp',root/'tests/CMakeLists.txt',root/'tests/crpropa/diagnose_network_photo_rates.py',root/'tests/crpropa/run_network_rate_scan.py']
    manifest=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),histories=0,seeds='none; deterministic means',host=platform.platform(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        source_status=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=root,text=True).strip(),
        tool_sha256={str(f.relative_to(root)):sha(f) for f in source},dependencies=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__),
        settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        input_sha256={str(f):sha(f) for f in [data/'manifest.json',data/'decay_manifest.json',data/'photon_manifest.json',reference/'shared_manifest.json']})
    def save(): (out/'run-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    try:
        pd=a.pd_photons=='On';photons=a.decay_photons=='On'
        verify(reference,'shared_manifest.json');raw=RawBundle(data,pd);decay=DecayBundle(data,reference,raw,a.provider,photons)
        bounds=energy_bounds(raw,decay);grid=grids(bounds,a.energy_per_nucleon,a.nodes)
        kernels={}
        for isotope,channels in decay.channels.items():
            if not channels:continue
            kernels[f'A{isotope[0]:02}Z{isotope[1]:02}']=[dict(code=c['code'],tau_s=c['tau'],
                products={f'A{p[0]:02}Z{p[1]:02}':n for p,n in c['kernel']['products'].items()},
                energy_per_nucleon_ratio=c['kernel']['ratio'],lepton_rewards_per_E_over_A=c['kernel']['leptons'],
                parent_boosted_photon_rest_lines_eV=c['photons']) for c in channels]
        (out/'decay-kernels.json').write_text(json.dumps(dict(kernels=kernels,electron_mass_eV=decay.me,mpc_seconds=decay.mpc_seconds,
            beta_sampler='51-node inverse CDF, uniform within each probability interval; isotropic mean boosts',
            free_neutron_decay=False,mass_MeV={f'A{a:02}Z{z:02}':decay.mass[a,z] for a,z in raw.reached if a>1}),indent=2,allow_nan=False)+'\n')
        grid_json={f'A{a:02}Z{z:02}':dict(bounds_E_over_A0=v,nodes_eV=grid[a,z].tolist()) for (a,z),v in bounds.items()}
        (out/'energy-grid.json').write_text(json.dumps(grid_json,indent=2)+'\n')
        slices=[];exports=[]
        if a.provider=='M':manifest['native_probe_sha256']=sha(a.native_probe.resolve());save()
        for j in range(max(map(len,grid.values()))):
            energies={p:float(values[min(j,len(values)-1)]) for p,values in grid.items()};native=None
            if a.provider=='M':
                gridfile=out/f'grid-{j:03}.tsv';gridfile.write_text(''.join(f'{a} {z} {energies[a,z]:.17g}\n' for a,z in raw.reached))
                native=out/f'native-{j:03}.tsv'
                command=[str(a.native_probe.resolve()),'--data',str(data),'--photon-table',str(reference/'shared_photons.tsv'),
                    '--energy-per-nucleon',str(a.energy_per_nucleon),'--output',str(native),'--photons',a.pd_photons,
                    '--background-step',str(a.background_step),'--eps-rel',str(a.eps_rel),'--combined','On','--decay-photons',a.decay_photons,'--grid',str(gridfile)]
                with (out/f'native-{j:03}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=log,check=True,timeout=600,env={**os.environ,'OMP_NUM_THREADS':'1'})
                exports.append(native_decay(native,photons));manifest.setdefault('native_commands',[]).append(command)
            model=build(data,reference,a.provider,a.energy_per_nucleon,pd,a.order,native,edges,bundle=raw,state_energies=energies)
            slices.append(model)
            (out/f'pd-slice-{j:03}.json').write_text(json.dumps(model,indent=2,allow_nan=False)+'\n')
            manifest['completed_coefficient_slices']=j+1;save()
        operator=assemble(slices,grid,decay,a.energy_per_nucleon,edges,pd,photons,exports if a.provider=='M' else None)
        save_npz(out/'operator.npz',operator['matrix'])
        result=solve(operator,a.distance_mpc,a.method,a.rtol,a.atol)
        (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        if a.path_control:
            other=solve(operator,a.distance_mpc,'Radau' if a.method=='BDF' else 'BDF',a.rtol/10,a.atol/10)
            (out/'path-control.json').write_text(json.dumps(other,indent=2,allow_nan=False)+'\n')
            manifest['path_reward_absolute_differences']={k:abs(v-other['rewards'][k]) for k,v in result['rewards'].items()}
        manifest.update(status='complete PD+decay deterministic diagnostic; no acceptance verdict',energy_states=len(operator['states']),finite_decay_rates=True,
            energy_remap='positive linear-energy number/mean preserving; DAG-derived bounds; no physical clipping',nuclear_photon_energy_debit=False)
        print(manifest['status'],manifest['energy_states'],result['initial_state_yield'],result['rewards']['photon_count'],result['checks'])
    except BaseException as error:manifest.update(status='failed or incomplete',error=str(error));raise
    finally:
        manifest['seconds']=time.monotonic()-start;manifest['outputs_sha256']={f.name:sha(f) for f in out.iterdir() if f.is_file() and f.name!='run-manifest.json'};save()


if __name__=='__main__':main()
