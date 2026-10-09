#!/usr/bin/env python3
"""Prepare a pinned primary photodisintegration network or tabulated masses.

All energy nodes and cross-section values are retained (MeV and mb). Outputs
belong outside Git. NIST/NuDat3 masses follow the pinned CRPropa generator,
without its approximate-mass fallback. No third-party Python modules are used.
"""
import argparse
from collections import defaultdict, deque
import hashlib
import json
import math
from pathlib import Path
import subprocess

REVISION = '34b10e27375c0e9894d7c984ac15fb96e12975ef'
TERMINALS = {(5, 2), (5, 3), (6, 2), (8, 4), (9, 5), (12, 5), (16, 9),
             (20, 11), (39, 20), (42, 21), (46, 23), (47, 24), (50, 25)}
EMITTED = ((1, 0), (1, 1), (2, 1), (3, 1), (3, 2), (4, 2))


def prepare(repo, output, policy, primary=(56, 26), mass_only=False):
    if primary[0] < 2 or not 0 < primary[1] <= primary[0]:
        raise ValueError("Invalid primary A,Z")
    def raw(name):
        return subprocess.check_output(['git', '-C', str(repo), 'show', f'{REVISION}:{name}'])
    head = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
    if head != REVISION:
        raise ValueError(f'Expected checkout {REVISION}, found {head}')
    sources = {}
    def text(name):
        data = raw(name)
        sources[name] = hashlib.sha256(data).hexdigest()
        return data.decode()
    def rows(name):
        return [line.split() for line in text(name).splitlines() if line.strip() and not line.startswith('#')]
    totals, channels, grids = {}, defaultdict(list), {}
    for label, directory, total, exclusive in (
            ('external', 'PD_external', 'xs_sum.txt', 'xs_excl.txt'),
            ('talys', 'PD_Talys1.8_Khan', 'xs_pd_sum.txt', 'xs_pd_thin.txt')):
        if mass_only: break
        base = 'tables/' + directory + '/'
        grid = [float(v) for row in rows(base + 'eps.txt') for v in row]
        if len(grid) < 2 or any(not math.isfinite(x) or x <= 0 for x in grid) or any(x >= y for x, y in zip(grid, grid[1:])):
            raise ValueError('Invalid photon-energy grid')
        grids[label] = grid
        for filename, is_channel in ((total, False), (exclusive, True)):
            for row in rows(base + filename):
                z, n = map(int, row[:2]); a = z + n
                if (a < 12) != (label == 'external'): continue
                offset = 3 if is_channel else 2
                values = [float(v) for v in row[offset:]]
                if len(values) != len(grid) or any(not math.isfinite(v) or v < 0 for v in values):
                    raise ValueError(f'Invalid cross section {a},{z}')
                if is_channel:
                    if not any(values): continue
                    code = row[2].zfill(6)
                    if len(code) != 6 or not code.isdigit(): raise ValueError('Invalid channel')
                    channels[a, z].append((code, values))
                else:
                    if (a, z) in totals: raise ValueError('Duplicate total')
                    totals[a, z] = values
    reached, pending = set(), deque([] if mass_only else [primary])
    while pending:
        a, z = key = pending.popleft()
        if key in reached: continue
        reached.add(key)
        if a == 1: continue
        if key in TERMINALS:
            if key in totals and channels[key]: raise ValueError('Terminal unexpectedly has complete data')
            continue
        if key not in totals or not channels[key]: raise ValueError(f'Missing cross sections for {key}')
        for code, _ in channels[key]:
            emitted = [(int(count), identity) for count, identity in zip(code, EMITTED) if int(count)]
            da = sum(count * identity[0] for count, identity in emitted)
            dz = sum(count * identity[1] for count, identity in emitted)
            residual = (a - da, z - dz)
            if da <= 0 or residual[0] < 1 or residual[1] < 0 or residual[1] > residual[0] or (residual[0] > 1 and residual[1] == 0):
                raise ValueError(f'Invalid residual {key} {code}: {residual}')
            pending.append(residual)
            pending.extend(identity for _, identity in emitted)
    if not mass_only and primary == (56,26) and TERMINALS - reached:
        raise ValueError('Terminal whitelist does not match reachable network')
    active_terminals = TERMINALS & reached
    # Reproduce constants and mass precedence from calc_mass.py/units.py.
    # Values are SI, converted to nuclear rest energy in MeV.
    amu, electron, c, mev_j = 1.660538921e-27, 9.10938291e-31, 2.99792458e8, 1.602176487e-13
    masses, provenance = {}, {}
    a = z = 0
    for line in text('tables/mass_NIST.txt').splitlines()[4:]:
        if line.startswith('Atomic Number'): z = int(line.split('=')[1])
        elif line.startswith('Mass Number'): a = int(line.split('=')[1])
        elif line.startswith('Relative Atomic Mass') and a > 1:
            value = line.split('=')[1].split('(')[0].strip()
            if value:
                masses[a, z] = (float(value) * amu - z * electron) * c*c / mev_j
                provenance[a, z] = 'NIST'
    for line in text('tables/decay_NuDat3.txt').splitlines()[1:]:
        cols = line.split('\t')
        try: z, n, excess = int(cols[2]), int(cols[3]), float(cols[6])
        except (ValueError, IndexError): continue
        key = (z+n, z)
        if key not in masses:
            masses[key] = ((z+n) * amu - z * electron) * c*c / mev_j + excess * 1e-3
            provenance[key] = 'NuDat3'
    for key in reached:
        if key[0] > 1 and (key not in masses or masses[key] <= 0 or not math.isfinite(masses[key])):
            raise ValueError(f'Missing tabulated mass {key}')
    text('units.py'); text('calc_mass.py')
    if not mass_only:
        text('calc_photodisintegration.py'); text('interactionRate.py')
    output.mkdir(parents=True, exist_ok=False)
    def write(name, lines):
        (output/name).write_text('\n'.join(lines) + '\n')
    def format_values(values): return ' '.join(format(v, '.17g') for v in values)
    if mass_only:
        keys = sorted(key for key in masses if key[0] > 1 and 0 < key[1] <= key[0]
                      and math.isfinite(masses[key]) and masses[key] > 0)
        write('masses.tsv', ['# A Z nuclear_mass_MeV source'] +
              [f'{a} {z} {masses[a,z]:.17g} {provenance[a,z]}' for a,z in keys])
        write('provenance.txt', ['mcray-talys-v1', REVISION, 'grid_policy=' + policy,
              'mass-only=true; photodisintegration=unavailable; nuclear-decays=unavailable',
              'mass=NIST>NuDat3; no-approximate-fallback; free-p,n=mcray'])
        manifest = {'revision': REVISION, 'source_sha256': sources, 'mass_only': True,
                    'mass_isotopes': len(keys),
                    'payload_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                      for p in sorted(output.iterdir())}}
        (output/'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        verify(output)
        print(f'Prepared {len(keys)} tabulated masses in {output}')
        return
    padded_grids = {}
    for label, grid in grids.items():
        step = (math.log(grid[-1]) - math.log(grid[0])) / (len(grid)-1)
        padded_grids[label] = grid + [grid[-1]*math.exp(step*j) for j in range(1, 514-len(grid))]
        write('grid_' + label + '.tsv', ['# K/MeV'] + [format(v, '.17g') for v in padded_grids[label]])
    def padded(a, values):
        return values + [0.] * (len(padded_grids['external' if a < 12 else 'talys'])-len(values))
    write('masses.tsv', ['# A Z nuclear_mass_MeV source'] + [f'{a} {z} {masses[a,z]:.17g} {provenance[a,z]}' for a,z in sorted(reached) if a > 1])
    write('terminals.tsv', ['# A Z: explicitly approved photodisintegration terminal states; nuclear decays disabled'] + [f'{a} {z}' for a,z in sorted(active_terminals)])
    write('total.tsv', ['# A Z sigma/mb at grid nodes'] + [f'{a} {z} ' + format_values(padded(a, totals[a,z])) for a,z in sorted(reached) if a > 1 and (a,z) not in TERMINALS])
    write('channels.tsv', ['# A Z channel(n,p,d,t,He3,He4) sigma/mb at grid nodes'] + [f'{a} {z} {code} ' + format_values(padded(a, values)) for a,z in sorted(reached) if a > 1 and (a,z) not in TERMINALS for code,values in channels[a,z]])
    write('provenance.txt', ['mcray-talys-v1', REVISION, 'grid_policy=' + policy,
          f'primary=A{primary[0]:02}Z{primary[1]:02}; nuclear_decays=Off; photon_emission=Off',
          'sources=PD_Talys1.8_Khan(full-total,thin-channels),PD_external(A<12)',
          'interpolation=log-log; zero-interval=sigma-versus-logK; outside-grid=zero',
          'mass=NIST>NuDat3; no-approximate-fallback; free-p,n=mcray',
          'terminals=' + ','.join(f'A{a:02}Z{z:02}' for a,z in sorted(active_terminals))])
    manifest = {'revision': REVISION, 'source_sha256': sources, 'reachable_states': sorted(reached), 'primary': list(primary),
                'payload_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir())}}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    verify(output)
    print(f'Prepared {len(reached)} states in {output}')


def verify(output):
    manifest = json.loads((output/'manifest.json').read_text())
    if manifest['revision'] != REVISION: raise ValueError('Wrong data revision')
    for name, expected in manifest['payload_sha256'].items():
        if hashlib.sha256((output/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Checksum mismatch: ' + name)
    print('All payload SHA-256 checksums match')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--grid-policy', choices=('crpropa-padding',), required=False)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--A', type=int, default=56, help='Primary mass number (default 56)')
    parser.add_argument('--Z', type=int, default=26, help='Primary charge number (default 26)')
    parser.add_argument('--mass-only', action='store_true', help='Export all available tabulated masses without cross sections')
    args = parser.parse_args()
    if args.verify: verify(args.output)
    elif not args.repo or not args.grid_policy: parser.error('--repo and --grid-policy are required for preparation')
    else: prepare(args.repo, args.output, args.grid_policy, (args.A,args.Z), args.mass_only)
