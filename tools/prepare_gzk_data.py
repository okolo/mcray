#!/usr/bin/env python3
"""Prepare a switchable Strict/CRPropa nuclear GZK bundle from pinned sources.

Requires an existing verified corrected decay bundle. No approximate PD cross
sections are introduced. Reachability and selected-policy checks run in mcray.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from prepare_talys_data import REVISION, TERMINALS, prepare as prepare_masses, verify
from prepare_nuclear_decay_data import verify as verify_decay


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(repo, decay, output):
    verify_decay(decay)
    if output.exists():
        raise ValueError('Output already exists; preserve previous bundles')
    if subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip() != REVISION:
        raise ValueError('Wrong source revision')
    sources = {}

    def raw(name):
        data = subprocess.check_output(['git', '-C', str(repo), 'show', f'{REVISION}:{name}'])
        sources[name] = hashlib.sha256(data).hexdigest()
        return [x.split() for x in data.decode().splitlines() if x.strip() and not x.startswith('#')]

    with tempfile.TemporaryDirectory(prefix='mcray-gzk-masses-') as temporary:
        mass_dir = Path(temporary) / 'masses'
        prepare_masses(repo, mass_dir, 'crpropa-padding', mass_only=True)
        measured = {}
        for line in (mass_dir/'masses.tsv').read_text().splitlines():
            if line and not line.startswith('#'):
                a, z, mass, source = line.split()
                if int(a) <= 56:
                    measured[int(a), int(z)] = (float(mass), source)
        sources.update(json.loads((mass_dir/'manifest.json').read_text())['source_sha256'])

    inventory_path = decay/'decay_sources/corrected/inventory.json'
    inventory = json.loads(inventory_path.read_text())
    masses, active, stable, fallback, unsupported = [], [], [], [], []
    for a in range(2, 57):
        for z in range(1, a+1):
            # Same pinned constants as calc_mass.py and prepare_talys_data.py.
            estimate = (a*1.660538921e-27-z*9.10938291e-31)*2.99792458e8**2/1.602176487e-13
            mass, source = measured.get((a,z), (estimate, 'CRPropaApprox'))
            masses.append(f'{a} {z} {mass:.17g} {source}')
            entries = inventory.get(f'{a},{z}')
            if not entries:
                raise ValueError(f'Missing source inventory for {a},{z}')
            if all(d['tau'] is None for d in entries):
                stable.append(f'{a} {z}')
            for d in entries:
                if d['tau'] is None:
                    continue
                if d['code'] is None:
                    # The generator also uses br=0 for unspecified single
                    # branches with valid lifetimes. Never drop those known
                    # modes. Only unsupported zero-weight entries are absent
                    # from the reference writer and irrelevant to its table.
                    if d['branch'] == 0:
                        continue
                    unsupported.append({'A':a, 'Z':z, 'mode':d['mode']})
                    continue
                row = f"{a} {z} {d['code'].zfill(5)} {d['tau']:.17g}"
                (fallback if d['tau'] < 1e-90 else active).append(row)

    payload = {}
    payload['masses.tsv'] = ['# A Z mass_MeV source; CRPropaApprox is opt-in'] + masses
    payload['decays.tsv'] = ['# A Z code partial_mean_lifetime_s; corrected known channels'] + active
    payload['fallback_decays.tsv'] = ['# Explicit CRPropa artificial dripping, enabled only by CRPropa policy'] + fallback
    payload['decay_stable.tsv'] = ['# Explicit inactive/stable states from corrected generator'] + stable
    payload['unsupported_decays.tsv'] = ['# A Z source-mode; calculation must reject reachable entries'] + [
        f"{d['A']} {d['Z']} {d['mode']}" for d in unsupported]
    payload['terminals.tsv'] = ['# Original explicitly approved PD terminals only'] + [f'{a} {z}' for a,z in sorted(TERMINALS)]
    totals, channels = [], []
    for label, directory, total, exclusive in (
            ('external','PD_external','xs_sum.txt','xs_excl.txt'),
            ('talys','PD_Talys1.8_Khan','xs_pd_sum.txt','xs_pd_thin.txt')):
        grid = [float(x) for row in raw(f'tables/{directory}/eps.txt') for x in row]
        import math
        step = math.log(grid[-1]/grid[0])/(len(grid)-1)
        padded = grid + [grid[-1]*math.exp(step*j) for j in range(1,514-len(grid))]
        payload[f'grid_{label}.tsv'] = ['# K/MeV'] + [f'{x:.17g}' for x in padded]
        for filename, is_channel in [(total,False),(exclusive,True)]:
            for row in raw(f'tables/{directory}/{filename}'):
                z,n = map(int,row[:2]); a=z+n
                if not 2<=a<=56 or not 0<z<=a or (a<12)!=(label=='external'):
                    continue
                if (a,z) not in measured:
                    raise ValueError(f'PD source lacks a verified mass: {a},{z}')
                offset=3 if is_channel else 2
                values=list(map(float,row[offset:]))
                if len(values)!=len(grid) or any(not math.isfinite(x) or x<0 for x in values):
                    raise ValueError('Invalid source cross-section grid')
                if is_channel and not any(values):
                    continue
                prefix=f'{a} {z}' + (f' {row[2].zfill(6)}' if is_channel else '')
                (channels if is_channel else totals).append(prefix+' '+
                    ' '.join(f'{x:.17g}' for x in values+[0.]*(513-len(values))))
    payload['total.tsv']=['# A Z full-total sigma/mb']+totals
    payload['channels.tsv']=['# A Z exclusive channel sigma/mb']+channels
    payload['provenance.txt']=['mcray-talys-v1',REVISION,'grid_policy=crpropa-padding',
        'gzk_compatibility=true','scope=A2..56; source-inventory=all; runtime-reachability=required',
        'mass=NIST>NuDat3>CRPropaApprox(policy-selected)',
        'PD=raw-pinned-full-total-thin-channels; missing-PD=error; original-terminals-only']
    payload['decay_provenance.txt']=['mcray-decay-v1',REVISION,
        'branching=exclusive-normalization-corrected; Be11-B-P=10010',
        'artificial-dripping=separate-fallback-table; activation=CRPropa-policy-only',
        'pure-p/n-spectators=prompt-free-nucleons; energy=by-A',
        'unsupported-source-modes=explicit-inventory; missing-decay=error']
    output.mkdir(parents=True)
    for name, lines in payload.items():
        (output/name).write_text('\n'.join(lines)+'\n')
    shutil.copyfile(inventory_path,output/'gzk_inventory.json')
    (output/'unsupported_decays.json').write_text(json.dumps(unsupported,indent=2)+'\n')
    decay_manifest=json.loads((decay/'decay_manifest.json').read_text())
    manifest={'revision':REVISION,'source_sha256':sources,
        'corrected_decay_source_sha256':decay_manifest['source_sha256'],
        'corrected_inventory_sha256':sha(inventory_path),'preparation_sha256':sha(Path(__file__)),
        'approximate_masses':sum(x.endswith('CRPropaApprox') for x in masses),
        'fallback_channels':len(fallback),'known_channels':len(active),'unsupported_modes':unsupported,
        'payload_sha256':{p.name:sha(p) for p in sorted(output.iterdir())}}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (output/'decay_manifest.json').write_text(json.dumps({'revision':REVISION,
        'source_sha256':decay_manifest['source_sha256'],
        'payload_sha256':{name:sha(output/name) for name in
            ['decays.tsv','decay_stable.tsv','fallback_decays.tsv','decay_provenance.txt','gzk_inventory.json','unsupported_decays.tsv']}},indent=2)+'\n')
    verify(output)
    print(f'Prepared switchable GZK bundle: {len(active)} known and {len(fallback)} fallback channels')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',required=True,type=Path)
    parser.add_argument('--decay-data',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    prepare(args.repo.resolve(),args.decay_data.resolve(),args.output.resolve())
