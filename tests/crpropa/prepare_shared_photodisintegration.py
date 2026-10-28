#!/usr/bin/env python3
"""Prepare CRPropa PD tables using its pinned integrator and shared photon density.

Only the audited Fe56 network is exported. Cross sections come from the verified
mcray raw bundle, with exactly the approved 513-node padding. NumPy/SciPy and
GitPython are preparation dependencies. Production physics code is unchanged.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from prepare_nuclear_decay_data import verify
from prepare_talys_data import REVISION
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo', type=Path, required=True)
parser.add_argument('--data', type=Path, required=True)
parser.add_argument('--photons', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
verify(args.data)
args.output.mkdir(parents=True, exist_ok=False)
source = args.output / 'sources'; source.mkdir()
sources = {}
for name in ('interactionRate.py', 'units.py', 'gitHelp.py', 'calc_photodisintegration.py'):
    raw = subprocess.check_output(['git', '-C', str(args.repo), 'show', f'{REVISION}:{name}'])
    (source / name).write_bytes(raw); sources[name] = hashlib.sha256(raw).hexdigest()
sys.path.insert(0, str(source.resolve()))
import interactionRate
from units import eV, Mpc
photons = np.loadtxt(args.photons)
energy, density = photons[:,0], photons[:,1]
if np.any(~np.isfinite(photons)) or np.any(energy<=0) or np.any(np.diff(energy)<=0) or np.any(density<0):
    raise ValueError('Invalid common photon table')
class SharedField:
    def getDensity(self, eps, z=0):
        x = np.asarray(eps)/eV
        j = np.clip(np.searchsorted(energy,x)-1,0,len(energy)-2)
        fraction = np.log(x/energy[j])/np.log(energy[j+1]/energy[j])
        positive = (density[j]>0)&(density[j+1]>0)
        value = np.where(positive, np.exp(np.log(density[j],where=density[j]>0,out=np.zeros_like(fraction))*(1-fraction)+
            np.log(density[j+1],where=density[j+1]>0,out=np.zeros_like(fraction))*fraction), density[j]*(1-fraction)+density[j+1]*fraction)
        return np.where((x>=energy[0])&(x<=energy[-1]),value,0.)*1e6/eV
field = SharedField()
gamma = np.logspace(6,14,201)
grids = {label: np.loadtxt(args.data/f'grid_{label}.tsv')*1e6*eV for label in ('talys','external')}
rate_rows, branch_rows, exclusive_rows = [], [], []
for file, channels in (('total.tsv',False),('channels.tsv',True)):
    grouped = defaultdict(list)
    for line in (args.data/file).read_text().splitlines():
        if not line or line.startswith('#'): continue
        row=line.split(); a,z=map(int,row[:2]); offset=3 if channels else 2
        xs=np.array(list(map(float,row[offset:])))*1e-31
        rate=interactionRate.calc_rate_eps(grids['external' if a<12 else 'talys'],xs,gamma,field)
        if np.any(~np.isfinite(rate)) or np.any(rate<0):raise ValueError('Invalid generated rate')
        if channels:
            grouped[a,z].append((row[2],rate));exclusive_rows.append((a,z,row[2],rate))
        else: rate_rows.append((z,a-z,rate))
    for (a,z),rows in grouped.items():
        total=np.sum([rate for _,rate in rows],axis=0)
        for code,rate in rows:
            branch_rows.append((z,a-z,code,np.divide(rate,total,out=np.zeros_like(rate),where=total>0)))
pd=args.output/'Photodisintegration';pd.mkdir()
with (pd/'rate_Shared.txt').open('w') as out:
    out.write('# Pinned CRPropa integrator; common photon table; gamma=10^6..10^14, 201 nodes\n')
    for z,n,rate in rate_rows:out.write(f'{z} {n} '+' '.join(f'{v:.9g}' for v in rate)+'\n')
with (pd/'branching_Shared.txt').open('w') as out:
    out.write('# Exclusive rates normalized by exclusive sum, as in the pinned generator\n')
    for z,n,code,branch in branch_rows:out.write(f'{z} {n} {code} '+' '.join(f'{v:.9g}' for v in branch)+'\n')
# The native constructor loads this file even when photon emission is disabled.
(pd/'photon_emission_Sha.txt').write_text('# Disabled diagnostic photon emission; no gamma channels requested\n')
with (args.output/'exclusive_rates.tsv').open('w') as out:
    out.write('# A Z channel rate_per_Mpc at gamma=10^6..10^14, 201 nodes\n')
    for a,z,code,rate in exclusive_rows:out.write(f'{a} {z} {code} '+' '.join(f'{v:.17g}' for v in rate)+'\n')
shutil.copyfile(args.photons,args.output/'shared_photons.tsv')
for name in ('nuclear_mass.txt','nuclear_decay.txt'):
    shutil.copyfile(args.data/'decay_sources/corrected/data'/name,args.output/name)
# Native TabularPhotonField stores E*dn/dE in m^-3, with energies in joules.
scaling=args.output/'Scaling';scaling.mkdir()
np.savetxt(scaling/'Shared_photonEnergy.txt',energy*eV,fmt='%.17g')
np.savetxt(scaling/'Shared_photonDensity.txt',energy*density*1e6,fmt='%.17g')
import scipy
manifest={'preparation_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'dependencies':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__},'revision':REVISION,'sources_sha256':sources,'photon_sha256':hashlib.sha256(args.photons.read_bytes()).hexdigest(),
 'decay_manifest_sha256':hashlib.sha256((args.data/'decay_manifest.json').read_bytes()).hexdigest(),
 'field':'Shared, static z=0; common log-log dn/dE interpolation; zero intervals linear in logE',
 'runtime_field_interpolation':'CRPropa TabularPhotonField is linear; PD transport uses preconvolved rates',
 'gamma_nodes':201,'payload_sha256':{str(p.relative_to(args.output)):hashlib.sha256(p.read_bytes()).hexdigest()
   for p in sorted(args.output.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}}
(args.output/'shared_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Prepared {len(rate_rows)} total and {len(branch_rows)} exclusive curves in {args.output}')
