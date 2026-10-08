#!/usr/bin/env python3
"""Run frozen diagnostic decay comparisons; long profiles require external Linux.

No acceptance verdict is inferred from the diagnostic z scores. The prepared
bundle and built image are inputs; no source checkout or data is modified.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from prepare_nuclear_decay_data import verify


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def jobs(data, profile, family):
    masses = {tuple(map(int, p[:2])): float(p[2]) for line in (data / 'masses.tsv').read_text().splitlines()
              if line and not line.startswith('#') for p in [line.split()]}
    rates = Counter()
    for line in (data / 'decays.tsv').read_text().splitlines():
        if line and not line.startswith('#'):
            a, z, code, tau = line.split(); rates[int(a), int(z)] += 1 / float(tau)
    batches, size = {'smoke': (2, 1000), 'pilot': (20, 1000), 'extended': (20, 10000)}[profile]
    cases = []
    for a, z in ((12, 5), (20, 11), (11, 4)):
        for depth in (0.1, 1., 3.):
            length = depth * 2.99792458e8 * 1e9 / rates[a, z] / 3.085677581491367e22
            cases.append((f'A{a:02}Z{z:02}-tau{depth:g}', a, z, length))
    for a, z in ((56, 26), (8, 4), (5, 2), (5, 3)):
        cases.append((f'A{a:02}Z{z:02}-prompt-or-stable', a, z, 1.))
    if family == 'talys':
        cases=[]
        for per_nucleon in (1e18,10**18.6,1e19):
            for length in (1.,10.):
                for process in ('talys','talys-decay'):
                    cases.append((f'Fe56-E{per_nucleon:.6g}-L{length:g}-{process}',56,26,length,per_nucleon,process))
    else:
        cases=[(*case,None,'decay') for case in cases]
    result = []
    for label, a, z, length, per_nucleon, process in cases:
        for batch in range(batches):
            for code in ('mcray', 'crpropa'):
                seed = (2026101200 if code == 'mcray' else 2026201200) + 1000 * len(result) + batch
                result.append({'id': f'{label}-{code}-b{batch:02}', 'case': label, 'code': code,
                    'a': a, 'z': z, 'gamma': per_nucleon*a/(masses[a,z]*1e6) if per_nucleon else 1e9,
                    'per_nucleon_eV': per_nucleon or 1e9 * masses[a,z] * 1e6 / a, 'process':process,
                    'length_mpc': length, 'histories': size, 'seed': seed, 'batch': batch,
                    'step_fraction': 0.03})
    return result


def command(job):
    if job['code'] == 'mcray':
        result = ['/src/build/Release/tests/nuclear_scenario', '--process', job['process'], '--talys-data', '/data',
            '--a', str(job['a']), '--charge', str(job['z']), '--energy-per-nucleon', str(job['per_nucleon_eV']),
            '--distance-mpc', str(job['length_mpc']), '--primaries', str(job['histories']),
            '--seed', str(job['seed']), '--particle-output', '--output', '/results/' + job['id'] + '.tsv']
        if job['process'] != 'decay': result += ['--photon-table','/shared/shared_photons.tsv']
        return result
    return ['env', 'CRPROPA_DATA_PATH=' + ('/data/decay_sources/corrected/data' if job['process']=='decay' else '/shared'),
        '/reference-build/crpropa_decay_probe', '--mode', 'transport', '--process', job['process'], '--a', str(job['a']),
        '--charge', str(job['z']), '--gamma', str(job['gamma']), '--distance-mpc', str(job['length_mpc']),
        '--step-fraction', str(job['step_fraction']), '--primaries', str(job['histories']),
        '--seed', str(job['seed']), '--output', '/results/' + job['id'] + '.tsv']


def summarize(path, job):
    histories = {}
    primary = f"A{job['a']:02}Z{job['z']:02}"
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'): continue
        b, h, species, weight, energy = line.split(); weight, energy = float(weight), float(energy)
        if weight != 1 or not math.isfinite(energy) or energy < 0:
            raise ValueError('Invalid product in ' + str(path))
        record = histories.setdefault((b, h), Counter())
        if species.startswith('A'):
            a, z = map(int, species[1:].split('Z'))
        else:
            a, z = {'proton': (1, 1), 'neutron': (1, 0), 'electron': (0, -1),
                    'positron': (0, 1), 'neutrinoE': (0, 0), 'neutrinoAE': (0, 0)}[species]
        record['baryons'] += a; record['charge'] += z
        record['survival'] += species == primary
        record['yield:' + species] += weight
        record['energy:' + species] += energy
        record['energy:all'] += energy
    if len(histories) != job['histories']:
        raise ValueError('Incomplete primary ensemble: ' + str(path))
    for record in histories.values():
        if record['baryons'] != job['a'] or record['charge'] != job['z']:
            raise ValueError('Baryon/charge violation: ' + str(path))
    keys = set().union(*(record.keys() for record in histories.values()))
    return {key: sum(record[key] for record in histories.values()) / len(histories) for key in keys}


def report(output, plan):
    grouped = {}
    for job in plan['jobs']:
        marker = output / (job['id'] + '.done.json')
        if not marker.exists(): continue
        entry = json.loads(marker.read_text())
        if digest(output / (job['id'] + '.tsv')) != entry['output_sha256']:
            raise ValueError('Completed output checksum changed: ' + job['id'])
        grouped.setdefault(job['case'], {}).setdefault(job['code'], []).append(entry['observables'])
    rows = []
    for case, codes in grouped.items():
        if set(codes) != {'mcray', 'crpropa'}: continue
        keys = set().union(*(sample.keys() for values in codes.values() for sample in values))
        for key in sorted(keys):
            if key in ('baryons', 'charge'): continue
            values = {code: [entry.get(key, 0.) for entry in entries] for code, entries in codes.items()}
            means = {code: statistics.mean(v) for code, v in values.items()}
            variances = {code: statistics.variance(v) / len(v) if len(v) > 1 else None for code, v in values.items()}
            se = math.sqrt(sum(variances.values())) if all(v is not None for v in variances.values()) else None
            difference = means['mcray'] - means['crpropa']
            rows.append({'case': case, 'observable': key, 'means': means, 'batch_counts': {c: len(v) for c, v in values.items()},
                         'difference': difference, 'combined_standard_error': se,
                         'diagnostic_z': difference / se if se else None,
                         'status': 'diagnostic-only; no acceptance verdict'})
    (output / 'comparison-diagnostic.json').write_text(json.dumps(rows, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--family', choices=('decay','talys'), default='decay')
    parser.add_argument('--reference-data', type=Path)
    parser.add_argument('--profile', choices=('smoke', 'pilot', 'extended'), default='smoke')
    parser.add_argument('--image', default='mcray-crpropa-validation:local')
    parser.add_argument('--external-server', action='store_true')
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--max-local-seconds', type=int, default=1800)
    args = parser.parse_args()
    if args.external_server and (platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'amd64')):
        parser.error('--external-server requires a native Linux amd64 host')
    if args.profile == 'extended' and not args.external_server and not args.plan_only:
        parser.error('Extended ensembles must run with the external Linux server script')
    if not 1 <= args.max_local_seconds <= 1800:
        parser.error('Local wall-clock limit must be between 1 and 1800 seconds')
    data, output = args.data.resolve(), args.output.resolve()
    if output == ROOT or ROOT in output.parents:
        parser.error('Results must be outside the Git repository')
    verify(data)
    shared = args.reference_data.resolve() if args.reference_data else None
    if args.family=='talys' and shared is None: parser.error('TALYS requires --reference-data with shared photon/rate tables')
    if shared:
        manifest=json.loads((shared/'shared_manifest.json').read_text())
        for name,checksum in manifest['payload_sha256'].items():
            if digest(shared/name)!=checksum: raise ValueError('Shared-data checksum mismatch: '+name)
    data_hashes = {name: digest(data / name) for name in ('manifest.json', 'decay_manifest.json')}
    if shared: data_hashes['shared_manifest.json']=digest(shared/'shared_manifest.json')
    image = json.loads(subprocess.check_output(['docker', 'image', 'inspect', args.image], text=True))[0]
    if image['Architecture'] != 'amd64' or image['Os'] != 'linux':
        parser.error('Reference image must be Linux amd64')
    image_id = image['Id']
    if args.resume:
        plan = json.loads((output / 'jobs.json').read_text())
        if plan['runner_sha256'] != digest(Path(__file__)):
            raise ValueError('Resume runner differs from the frozen plan')
        if plan['data_manifest_sha256'] != data_hashes or plan['profile'] != args.profile or plan.get('family','decay') != args.family or plan['image_id'] != image_id:
            raise ValueError('Resume inputs differ from the frozen plan')
    else:
        output.mkdir(parents=True, exist_ok=False)
        plan = {'profile': args.profile, 'family': args.family, 'data_manifest_sha256': data_hashes, 'image_id': image_id,
                'runner_sha256': digest(Path(__file__)), 'jobs': jobs(data, args.profile, args.family),
                'status': 'diagnostic ensemble; acceptance family not frozen'}
        (output / 'jobs.json').write_text(json.dumps(plan, indent=2) + '\n')
    if args.plan_only:
        print(f"Prepared {len(plan['jobs'])} jobs in {output}; no calculation started")
        return
    start = time.monotonic()
    for index, job in enumerate(plan['jobs']):
        if (output / (job['id'] + '.done.json')).exists(): continue
        elapsed = time.monotonic() - start
        if not args.external_server and elapsed >= args.max_local_seconds:
            report(output, plan)
            raise RuntimeError('Local wall limit reached; resume remaining jobs on the external Linux server')
        container_name = 'mcray-decay-' + hashlib.sha256(str(output).encode()).hexdigest()[:12] + '-' + str(index)
        cmd = ['docker', 'run', '--rm', '--name', container_name, '--platform', 'linux/amd64',
               '--mount', f'type=bind,source={data},target=/data,readonly',
               '--mount', f'type=bind,source={output},target=/results'] + ([] if shared is None else ['--mount',f'type=bind,source={shared},target=/shared,readonly']) + [image_id] + command(job)
        print(f"[{index+1}/{len(plan['jobs'])}] {job['id']}", flush=True)
        with (output / (job['id'] + '.log')).open('w') as log:
            try:
                subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT,
                               timeout=None if args.external_server else max(1., args.max_local_seconds-elapsed))
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                subprocess.run(['docker', 'stop', '--time', '2', container_name], stdout=log, stderr=subprocess.STDOUT)
                report(output, plan)
                raise
        observables = summarize(output / (job['id'] + '.tsv'), job)
        (output / (job['id'] + '.done.json')).write_text(json.dumps({
            'observables': observables, 'output_sha256': digest(output / (job['id'] + '.tsv')),
            'command': cmd, 'elapsed_seconds_since_start': time.monotonic() - start}, indent=2) + '\n')
    report(output, plan)
    print(f'Completed in {time.monotonic()-start:.1f} seconds; diagnostic only', flush=True)


if __name__ == '__main__':
    main()
