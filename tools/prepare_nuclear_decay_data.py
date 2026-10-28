#!/usr/bin/env python3
"""Extend a verified TALYS bundle with pinned, corrected CRPropa decay data.

Run in an isolated Python environment with NumPy, SciPy and GitPython. Generated
sources, original reference tables, corrected tables and hashes stay outside Git.
Only the approved branching normalization and Be11 B-P mapping are changed.
Unknown decay modes, artificial drip lifetimes, missing states and masses fail.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from prepare_talys_data import REVISION, verify as verify_talys

FILES = ('calc_mass.py', 'calc_decay.py', 'units.py', 'gitHelp.py',
         'tables/mass_NIST.txt', 'tables/decay_NuDat3.txt', 'tables/gamma_NuDat3.txt')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corrected_script(source):
    old = '            # normalize all branching ratios\n            for d in dList:\n                d.br /= brSum'
    new = '            # Recompute after converting inclusive branches to exclusive.\n            brSum = sum(d.br for d in dList)\n            for d in dList:\n                d.br /= brSum'
    if source.count(old) != 1 or source.count("'B-N': '10001',") != 1:
        raise ValueError('Pinned decay generator patch context changed')
    return source.replace(old, new).replace("'B-N': '10001',", "'B-N': '10001',\n        'B-P': '10010',")


def verify(directory):
    verify_talys(directory)
    manifest = json.loads((directory / 'decay_manifest.json').read_text())
    if manifest['revision'] != REVISION:
        raise ValueError('Wrong decay source revision')
    for name, digest in manifest['payload_sha256'].items():
        if sha(directory / name) != digest:
            raise ValueError('Decay checksum mismatch: ' + name)
    print('All nuclear decay SHA-256 checksums match')


def prepare(repo, talys, output):
    verify_talys(talys)
    if subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip() != REVISION:
        raise ValueError('Wrong source checkout revision')
    if output.exists():
        raise ValueError('Output directory already exists; keep previous data immutable')
    shutil.copytree(talys, output)
    reference = output / 'decay_sources' / 'reference'
    corrected = output / 'decay_sources' / 'corrected'
    sources = {}
    for name in FILES:
        data = subprocess.check_output(['git', '-C', str(repo), 'show', f'{REVISION}:{name}'])
        sources[name] = hashlib.sha256(data).hexdigest()
        for folder in (reference, corrected):
            path = folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    path = corrected / 'calc_decay.py'
    path.write_text(corrected_script(path.read_text()))
    # Child processes isolate the generator's globals and module import paths.
    # A separate inventory distinguishes explicitly stable isotopes from omissions.
    runner = '''import contextlib, json, runpy
with open('generator.log','w') as log, contextlib.redirect_stdout(log):
    runpy.run_path('calc_mass.py')
    env=runpy.run_path('calc_decay.py')
inventory={}
for z, row in enumerate(env['decayTable']):
    for n, channels in enumerate(row):
        if z+n>56 or z==0: continue
        inventory[f'{z+n},{z}']=[{'mode':d.mode,'tau':None if d.isStable() else float(d.tau),
            'code':env['modeDict'].get(d.mode), 'branch':float(d.br)} for d in channels]
with open('inventory.json','w') as out: json.dump(inventory,out,indent=2)
'''
    for folder in (reference, corrected):
        (folder / 'run.py').write_text(runner)
        subprocess.run([sys.executable, 'run.py'], cwd=folder, check=True,
                       env=dict(os.environ, GIT_DIR=str(repo / '.git')))
    inventory = json.loads((corrected / 'inventory.json').read_text())
    masses = {tuple(map(int, line.split()[:2])) for line in (output / 'masses.tsv').read_text().splitlines()
              if line.strip() and not line.startswith('#')}
    decays, stable = [], []
    for a, z in sorted(masses):
        entries = inventory.get(f'{a},{z}')
        if not entries:
            raise ValueError(f'Missing decay information: {a},{z}')
        if len(entries) > 1 and abs(sum(d['branch'] for d in entries) - 1.) > 1e-10:
            raise ValueError(f'Exclusive branching fractions do not sum to one: {a},{z}')
        active = [d for d in entries if d['tau'] is not None]
        if not active:
            stable.append(f'{a} {z}')
        for d in active:
            if d['tau'] <= 1e-90 or d['code'] is None:
                raise ValueError(f'Unapproved missing-data fallback or unsupported mode: {a},{z},{d}')
            decays.append(f"{a} {z} {d['code'].zfill(5)} {d['tau']:.17g}")
    (output / 'decays.tsv').write_text('# A Z channel(beta-,beta+,alpha,p,n) partial_mean_lifetime_s\n' + '\n'.join(decays) + '\n')
    (output / 'decay_stable.tsv').write_text('# A Z: stable or beta+ forbidden for fully stripped ions\n' + '\n'.join(stable) + '\n')
    (output / 'decay_provenance.txt').write_text('\n'.join([
        'mcray-decay-v1', REVISION,
        'branching=exclusive-normalization-corrected; Be11-B-P=10010',
        'ionization=CRPropa-electron-capture-correction; forbidden-beta-plus=stable',
        'kinematics=CRPropa; beta=same-gamma-no-recoil; emission=energy-by-A',
        'electrons=On; neutrinos=On; decay-photons=Off',
        'artificial-dripping=forbidden; missing-decay-data=error',
        'reference=decay_sources/reference/data/nuclear_decay.txt',
        'TALYS-provenance-describes-photodisintegration-only; decay-activation=explicit']) + '\n')
    manifest = {'revision': REVISION, 'source_sha256': sources, 'preparation_sha256': sha(Path(__file__)),
                'python': sys.version, 'dependencies': {name: version(name) for name in ('numpy', 'scipy', 'GitPython')},
                'decay_channels': len(decays), 'stable_states': len(stable),
                'payload_sha256': {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*'))
                    if p.is_file() and ('decay_sources' in p.parts or p.name.startswith('decay'))
                    and '__pycache__' not in p.parts}}
    (output / 'decay_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    verify(output)
    print(f'Prepared {len(decays)} decay channels and {len(stable)} stable states in {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--talys-data', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify(args.output.resolve())
    elif not args.repo or not args.talys_data:
        parser.error('--repo and --talys-data are required')
    else:
        prepare(args.repo.resolve(), args.talys_data.resolve(), args.output.resolve())
