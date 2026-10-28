#!/usr/bin/env python3
"""Add pinned raw PD photon cross sections and corrected decay lines to a bundle.

Numerical payloads stay outside Git. Missing transition coverage is retained
as a runtime/preflight error; source inconsistencies are never normalized.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from prepare_talys_data import REVISION, verify as verify_talys
from prepare_nuclear_decay_data import verify as verify_decay


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(repo, bundle, output, decay_reference):
    verify_talys(bundle)
    if output.exists():
        raise ValueError('Output exists; preserve immutable bundles')
    if subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()!=REVISION:
        raise ValueError('Wrong source revision')
    sources={}
    def raw(name):
        data=subprocess.check_output(['git','-C',str(repo),'show',f'{REVISION}:{name}'])
        sources[name]=hashlib.sha256(data).hexdigest()
        return [s.split() for s in data.decode().splitlines() if s.strip() and not s.startswith('#')]
    payload={}
    exceptions=[(2,1,1,1),(3,1,1,1),(3,1,2,1),(3,2,1,1),(3,2,2,1),(4,2,2,1),(4,2,3,2),(7,3,6,2),(8,3,7,3),(8,5,7,4)]
    payload["pd_photon_exceptions.tsv"]="# Explicit approved no-photon transitions; missing source record is not physical zero\n"+"\n".join(" ".join(map(str,key)) for key in exceptions)+"\n"
    grid=[float(x) for row in raw('tables/PD_Talys1.8_Khan/eps.txt') for x in row]
    for filename,target,offset in [('xs_photon_sum.txt','pd_photon_transitions.tsv',4),
                                    ('xs_photon_thin.txt','pd_photons.tsv',5)]:
        rows=[]
        for row in raw('tables/PD_Talys1.8_Khan/'+filename):
            z,n,zd,nd=map(int,row[:4]); a=z+n; ar=zd+nd
            if not 2<=a<=56: continue
            if not 1<=ar<a or not 0<=zd<=ar: raise ValueError('Invalid photon transition')
            values=list(map(float,row[offset:]))
            if len(values)!=len(grid) or len(values)>513 or any(not math.isfinite(x) or x<0 for x in values):
                raise ValueError('Invalid photon cross sections')
            prefix=f'{a} {z} {ar} {zd}'
            if offset==5:
                energy=float(row[4])
                if not math.isfinite(energy) or energy<=0: raise ValueError('Invalid line energy')
                prefix+=f' {energy:.17g}'
            rows.append(prefix+' '+' '.join(f'{x:.17g}' for x in values+[0.]*(513-len(values))))
        payload[target]='# A Z residual_A residual_Z '+('energy_MeV ' if offset==5 else '')+'sigma_mb[513]\n'+'\n'.join(rows)+'\n'
    # Hash the generator and its thinning rules as well as the raw payloads.
    for name in ['calc_photodisintegration.py','tables/PD_Talys1.8_Khan/A3_process_photonemission.py']:
        raw(name)
    reference_hashes={}
    if (bundle/'decays.tsv').exists():
        verify_decay(bundle)
        verify_decay(decay_reference)
        ref=decay_reference/'decay_sources/corrected/data/nuclear_decay.txt'
        reference_hashes['corrected_nuclear_decay.txt']=sha(ref)
        reference_hashes['decay_manifest.json']=sha(decay_reference/'decay_manifest.json')
        reference_hashes['gamma_NuDat3.txt']=sha(decay_reference/'decay_sources/corrected/tables/gamma_NuDat3.txt')
        entries={}
        for row in ref.read_text().splitlines():
            if not row.strip() or row.startswith('#'):continue
            s=row.split(); z,n=map(int,s[:2]); key=(z+n,z,s[2]); values=list(map(float,s[4:]))
            if len(values)%2:raise ValueError('Invalid corrected photon pairs')
            if key in entries and entries[key]!=values:raise ValueError('Ambiguous corrected decay photon channel')
            entries[key]=values
        coverage=[]; lines=[]
        for filename in ['decays.tsv','fallback_decays.tsv']:
            if not (bundle/filename).exists():continue
            for row in (bundle/filename).read_text().splitlines():
                if not row.strip() or row.startswith('#'):continue
                s=row.split();a,z=map(int,s[:2]);key=(a,z,s[2])
                if key not in entries:raise ValueError(f'Missing corrected photon coverage: {key}')
                coverage.append(f'{a} {z} {s[2]}')
                values=entries[key]
                for energy,prob in zip(values[::2],values[1::2]):
                    if not math.isfinite(energy) or energy<=0 or not math.isfinite(prob) or not 0<=prob<=1:
                        raise ValueError(f'Invalid photon probability: {key}, {energy}, {prob}')
                    lines.append(f'{a} {z} {s[2]} {energy:.17g} {prob:.17g}')
        payload['decay_photon_channels.tsv']='# Explicit reviewed channel coverage, including no-line channels\n'+'\n'.join(coverage)+'\n'
        payload['decay_photons.tsv']='# A Z channel energy_keV probability; attached once per selected channel\n'+'\n'.join(lines)+'\n'
    payload['photon_provenance.txt']='\n'.join(['mcray-nuclear-photons-v1',REVISION,
        'PD=native-convolution; denominator=pinned-xs_photon_sum; numerator=pinned-xs_photon_thin',
        'PD-missing-transition=error-except-ten-approved-light-transitions; see-pd_photon_exceptions.tsv',
        'PD-explicit-empty-transition=no-retained-lines-in-pinned-model; not-physical-zero',
        'decay=corrected-exclusive-channel; photons-once-before-compound-sequence',
        'boost=CRPropa-ultrarelativistic-isotropic; residual-and-beta-Q=not-debited',
        'line-correlation=independent-Bernoulli; exact-energy-conservation=not-claimed'])+'\n'
    shutil.copytree(bundle,output)
    for name,content in payload.items():(output/name).write_text(content)
    manifest={'revision':REVISION,'source_sha256':sources,'corrected_decay_reference_sha256':reference_hashes,
              'input_manifest_sha256':sha(bundle/'manifest.json'),'preparation_sha256':sha(Path(__file__)),
              'payload_sha256':{name:sha(output/name) for name in payload}}
    (output/'photon_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Prepared photon payload in {output}')


def verify(output):
    verify_talys(output)
    manifest=json.loads((output/'photon_manifest.json').read_text())
    if manifest['revision']!=REVISION:raise ValueError('Wrong photon source revision')
    for name,digest in manifest['payload_sha256'].items():
        if sha(output/name)!=digest:raise ValueError('Photon checksum mismatch: '+name)
    if (output/'decay_manifest.json').exists():verify_decay(output)
    print('All photon SHA-256 checksums match')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path);p.add_argument('--data',type=Path)
    p.add_argument('--decay-reference',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true');args=p.parse_args()
    if args.verify:verify(args.output.resolve())
    elif not args.repo or not args.data:p.error('--repo and --data required')
    else:
        prepare(args.repo.resolve(),args.data.resolve(),args.output.resolve(),
                (args.decay_reference or args.data).resolve())
        verify(args.output.resolve())
