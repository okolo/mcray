#!/usr/bin/env python3
"""Run the frozen initial-point acceptance once, with independent CPU jobs.

Insufficient precision is inconclusive. Never rerun samples until they pass.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
from run_transport_photon_pilots import sha, summarize, PIN, ROOT
from run_network_rate_scan import verify

PROTOCOL = Path(__file__).parent/'protocols/transport_photons_v1.json'
ENV = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}


def jobs(protocol):
    result = []
    for bg_index, bg in enumerate(protocol['backgrounds']):
        case = dict(id=bg, background=bg, process='talys-decay', pd_photons=True,
                    decay_photons=True, length_mpc=protocol['distance_mpc'],
                    energy_per_nucleon_eV=protocol['energy_per_nucleon_eV'])
        for code_index, code in enumerate(['mcray', 'crpropa']):
            for level in range(3):
                count = protocol['main_batches'] if level == 0 else protocol['refinement_batches']
                for batch in range(count):
                    seed = protocol['seed_base']+1000000*code_index+10000*bg_index+batch
                    result.append(dict(id=f'{bg}-{code}-l{level}-b{batch:02}', case=case,
                        code=code, level=level, batch=batch, histories=protocol['histories_per_batch'], seed=seed))
    total = (protocol['main_batches']+2*protocol['refinement_batches'])*protocol['histories_per_batch']
    if total+protocol['prior_pilot_histories']+protocol['other_controls_reservation'] > protocol['history_cap']:
        raise ValueError('History cap exceeded')
    return result


def summary(path, job, protocol):
    metadata = {}
    photon_modes = []
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            metadata.update(token.split('=', 1) for token in line.split() if '=' in token)
            if line.startswith('# PD_photons='):
                photon_modes.append(dict(token.split('=', 1) for token in line.split() if '=' in token))
    if (metadata.get('seed') != str(job['seed']) or metadata.get('process') != 'talys-decay' or
            float(metadata.get('L_mpc', 'nan')) != protocol['distance_mpc'] or
            not photon_modes or any(v.get('PD_photons') not in ['On', '1'] or
                v.get('decay_photons') not in ['On', '1'] for v in photon_modes)):
        raise ValueError('Output seed, scope or photon switches differ from frozen job')
    if job['code'] == 'mcray':
        control = protocol['native_controls'][job['level']]
        if (float(metadata.get('accuracy', 'nan')) != control['eps'] or
                float(metadata.get('background_step_decades', 'nan')) != control['step']):
            raise ValueError('Native numerical controls changed')
    elif float(metadata.get('step_fraction', 'nan')) != protocol['reference_step_fractions'][job['level']]:
        raise ValueError('Reference step control changed')
    values = summarize(path, job)
    extra = {'energy:bound': 0.}
    for i in range(4):
        extra[f'photon_bin_{i}_count'] = 0.
        extra[f'photon_bin_{i}_energy_eV'] = 0.
    edges = protocol['photon_bin_edges_eV']
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        b, h, species, weight, energy = line.split()
        energy = float(energy)
        if species.startswith('A'):
            extra['energy:bound'] += energy/job['histories']
        if species == 'photon':
            k = next(i for i in range(4) if edges[i+1] is None or energy < edges[i+1])
            extra[f'photon_bin_{k}_count'] += 1/job['histories']
            extra[f'photon_bin_{k}_energy_eV'] += energy/job['histories']
    values.update(extra)
    keys = protocol['targets'][job['case']['background']+'/'+job['code']]['means']
    return {key: values.get(key, 0.) for key in keys}


def interval(values, target, normalization, epsilon, spectral, protocol, family_size):
    from scipy.stats import t
    mean = statistics.mean(values)
    se = statistics.stdev(values)/math.sqrt(len(values))
    critical = float(t.ppf(1-(1-protocol['simultaneous_confidence'])/(2*family_size), len(values)-1))
    half = critical*se
    delta = normalization*(protocol['spectral_tolerance'] if spectral else protocol['integral_tolerance'])
    difference = mean-target
    status = 'inconclusive'
    # Constant/zero samples do not justify zero-width empirical intervals.
    if se > 0 and normalization > 0 and math.isfinite(half):
        if abs(difference)-half-epsilon > delta:
            status = 'disagreement'
        elif (abs(difference)+half+epsilon <= delta and
              (spectral or se <= .01*normalization)):
            status = 'agreement'
    return dict(mean=mean, target=target, normalization=normalization, difference=difference,
                standard_error=se, half_width=half, numerical_envelope=epsilon,
                tolerance=delta, status=status)


def report(out, protocol, all_jobs):
    grouped = {}
    for job in all_jobs:
        marker = out/(job['id']+'.done.json')
        if not marker.exists():
            continue
        done = json.loads(marker.read_text())
        if done['job'] != job or sha(out/(job['id']+'.tsv')) != done['output_sha256']:
            raise ValueError('Changed completed batch')
        actual = summary(out/(job['id']+'.tsv'), job, protocol)
        if actual != done['observables']:
            raise ValueError('Cached observables differ from complete raw histories')
        grouped.setdefault((job['case']['background'], job['code'], job['level']), []).append(actual)
    keys = list(protocol['targets']['cmb/mcray']['means'])
    # Include main and both predeclared refinement families in the union bound.
    family = len(keys)*len(protocol['backgrounds'])*2*3
    rows = []
    for (bg, code, level), batches in sorted(grouped.items()):
        expected_batches = protocol['main_batches'] if level == 0 else protocol['refinement_batches']
        if len(batches) != expected_batches:
            continue
        target = protocol['targets'][bg+'/'+code]
        normalization = protocol['targets'][bg+'/independent']['means']
        for key in keys:
            row = interval([v[key] for v in batches], target['means'][key], abs(normalization[key]),
                target['numerical_envelope'][key], key.startswith('photon_bin_'), protocol, family)
            row.update(background=bg, code=code, level=level, observable=key, batches=len(batches))
            rows.append(row)
    main = [r for r in rows if r['level'] == 0]
    complete = len(main) == family//3
    status = ('disagreement' if any(r['status'] == 'disagreement' for r in main) else
              'agreement' if complete and all(r['status'] == 'agreement' for r in main) else 'inconclusive')
    execution_complete = (len(grouped) == 12 and
        sum(len(v) for v in grouped.values()) == len(all_jobs))
    overall = ('disagreement' if any(r['status'] == 'disagreement' for r in rows) else
               'agreement' if execution_complete and all(r['status'] == 'agreement' for r in rows)
               else 'inconclusive')
    coefficients = []
    for bg in protocol['backgrounds']:
        reference = protocol['targets'][bg+'/independent']
        for code in ['mcray', 'crpropa']:
            native = protocol['targets'][bg+'/'+code]
            for key in keys:
                normalization = abs(reference['means'][key])
                delta = normalization*(protocol['spectral_tolerance'] if key.startswith('photon_bin_') else protocol['integral_tolerance'])
                epsilon = native['numerical_envelope'][key]+reference['numerical_envelope'][key]
                difference = native['means'][key]-reference['means'][key]
                coefficients.append(dict(background=bg, code=code, observable=key,
                    difference=difference, normalization=normalization, numerical_envelope=epsilon,
                    status=('agreement' if normalization and abs(difference)+epsilon <= delta else
                            'disagreement' if normalization and abs(difference)-epsilon > delta else 'inconclusive'),
                    role='native coefficient check' if code == 'mcray' else 'official representation diagnostic; not MC rejection'))
    native_coefficients = [r for r in coefficients if r['code'] == 'mcray']
    if any(r['status'] == 'disagreement' for r in native_coefficients):
        overall = 'disagreement'
    elif overall == 'agreement' and any(r['status'] != 'agreement' for r in native_coefficients):
        overall = 'inconclusive'
    result = dict(status=overall, execution_complete=execution_complete,
        deterministic_coefficient_comparisons=coefficients,
        scope=protocol['scope'], protocol_sha256=sha(PROTOCOL), simultaneous_family_size=family,
        main_family_status=status, observables=rows,
        note='Refinement rows are separate fixed controls, not empirical mean corrections. '
             'Numerical envelopes are observed deterministic grid/path differences, not rigorous global bounds. '
             'An inconclusive or partial report is never physical acceptance or a Release verdict.')
    (out/'acceptance-report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def command(job, protocol, bundle, references, native, reference_probe, target):
    common = ['--process', 'talys-decay', '--a', '56', '--charge', '26',
        '--distance-mpc', str(protocol['distance_mpc']), '--primaries', str(job['histories']),
        '--seed', str(job['seed']), '--output', str(target),
        '--photodisintegration-photons', 'On', '--decay-photons', 'On']
    bg = job['case']['background']
    if job['code'] == 'mcray':
        controls = protocol['native_controls'][job['level']]
        return ([str(native), '--talys-data', str(bundle/'nuclei'), '--photon-table',
            str(references/bg/'shared_photons.tsv'), '--particle-output', '--pd-photon-yields', 'CRPropa',
            '--energy-per-nucleon', str(protocol['energy_per_nucleon_eV']),
            '--background-step-decades', str(controls['step']), '--accuracy', str(controls['eps'])]+common, ENV)
    masses = {tuple(map(int, r[:2])): float(r[2]) for line in (bundle/'nuclei/masses.tsv').read_text().splitlines()
              if line and not line.startswith('#') for r in [line.split()]}
    gamma = 56*protocol['energy_per_nucleon_eV']/(masses[56, 26]*1e6)
    return ([str(reference_probe), '--mode', 'transport', '--gamma', str(gamma),
        '--step-fraction', str(protocol['reference_step_fractions'][job['level']])]+common,
        {**ENV, 'CRPROPA_DATA_PATH': str(references/bg)})


def check_coefficients(entry, protocol, cache, bundle, references, out, native_probe=None):
    case = cache/entry['case']
    for label in ['grid', 'native']:
        if sha(case/f'{label}-{entry["node"]:03}.tsv') != entry[label+'_sha256']:
            raise ValueError('Changed cached coefficient control')
    output = out/f'coefficient-{entry["background"]}-{entry["node"]:03}.tsv'
    controls = protocol['native_controls'][0]
    probe = native_probe or ROOT/'build/Release/tests/deterministic_pd_export'
    argv = [str(probe), '--data', str(bundle/'nuclei'),
        '--photon-table', str(references/entry['background']/'shared_photons.tsv'),
        '--energy-per-nucleon', str(protocol['energy_per_nucleon_eV']), '--output', str(output),
        '--photons', 'On', '--combined', 'On', '--decay-photons', 'On',
        '--background-step', str(controls['step']), '--eps-rel', str(controls['eps']),
        '--grid', str(case/f'grid-{entry["node"]:03}.tsv')]
    with output.with_suffix('.log').open('w') as log:
        subprocess.run(argv, check=True, env=ENV, stdout=log, stderr=log, timeout=600)
    def tokens(path):
        return [line.split() for line in path.read_text().splitlines() if line and not line.startswith('#')]
    old, new = tokens(case/f'native-{entry["node"]:03}.tsv'), tokens(output)
    if len(old) != len(new):
        raise ValueError('Coefficient inventory changed')
    for a, b in zip(old, new):
        if len(a) != len(b):
            raise ValueError('Coefficient row changed')
        for x, y in zip(a, b):
            try:
                equal = math.isclose(float(x), float(y), rel_tol=5e-12, abs_tol=0.)
            except ValueError:
                equal = x == y
            if not equal:
                raise ValueError('Current producer differs from frozen coefficients: '+str((a, b)))
    return dict(control=entry, command=argv, output_sha256=sha(output), status='matched')


def reserve_campaign(ledger, reservation):
    try:
        with ledger.open('x') as stream:
            json.dump(reservation, stream, indent=2); stream.write('\n')
    except FileExistsError:
        if json.loads(ledger.read_text()) != reservation:
            raise ValueError('Campaign already reserved at another output/source; audit history budget, do not rerun')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle', required=True, type=Path)
    p.add_argument('--references', required=True, type=Path)
    p.add_argument('--deterministic-cache', required=True, type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--expected-protocol-sha')
    p.add_argument('--campaign-ledger', type=Path)
    p.add_argument('--workers', type=int, default=min(8, os.cpu_count() or 1))
    p.add_argument('--external-server', action='store_true')
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    if a.workers < 1:
        p.error('--workers must be positive')
    if not a.plan_only and a.campaign_ledger is None:
        p.error('Execution requires --campaign-ledger to prevent duplicate campaigns')
    if not a.plan_only and (not a.external_server or platform.system() != 'Linux' or platform.machine() != 'x86_64'):
        p.error('Execution requires an external native Linux amd64 server')
    if not a.plan_only and a.expected_protocol_sha != sha(PROTOCOL):
        p.error('Launcher protocol hash is missing or differs from fetched source')
    protocol = json.loads(PROTOCOL.read_text()); all_jobs = jobs(protocol)
    bundle, references, cache = a.bundle.resolve(), a.references.resolve(), a.deterministic_cache.resolve()
    inputs = json.loads((bundle/'inputs-manifest.json').read_text())
    if inputs['data_commit'] != '34b10e27375c0e9894d7c984ac15fb96e12975ef' or inputs['reference_code'] != PIN:
        raise ValueError('Transferred input versions changed')
    for name, digest in inputs['payload_sha256'].items():
        if sha(bundle/name) != digest:
            raise ValueError('Transferred input checksum mismatch: '+name)
    for name, digest in protocol['data_manifest_sha256'].items():
        if sha(bundle/'nuclei'/name) != digest:
            raise ValueError('Frozen nuclear data changed')
    for bg, digest in protocol['prepared_reference_manifest_sha256'].items():
        if sha(references/bg/'shared_manifest.json') != digest:
            raise ValueError('Prepared control differs from frozen reference')
        verify(references/bg, 'shared_manifest.json')
    for name, digest in protocol['reference_files_sha256'].items():
        if sha(cache/name) != digest:
            raise ValueError('Frozen deterministic result changed: '+name)
    probe = references/'reference-build/crpropa_decay_probe'
    if sha(probe) != protocol['crpropa_probe_sha256']:
        raise ValueError('Pinned CRPropa binary changed')
    build = json.loads((references/'reference-build.json').read_text())
    if build['commit'] != PIN or build['probe_sha256'] != sha(probe):
        raise ValueError('Pinned CRPropa source provenance mismatch')
    native = ROOT/'build/Release/tests/nuclear_scenario'
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    status = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT, text=True).strip()
    if status and not a.plan_only:
        raise ValueError('Acceptance requires a clean frozen source checkout')
    parent = Path(os.environ['MCRAY_RESULTS_DIR']) if 'MCRAY_RESULTS_DIR' in os.environ else None
    out = a.output or (parent/'transport-photon-acceptance' if parent else None)
    if out is None:
        p.error('Supply --output or MCRAY_RESULTS_DIR')
    out = out.resolve()
    identity = dict(protocol_sha256=sha(PROTOCOL), source_commit=source,
        source_patch_sha256=__import__('hashlib').sha256(subprocess.check_output(['git','diff','HEAD'],cwd=ROOT)).hexdigest(),
        driver_sha256=sha(Path(__file__)), native_sha256=sha(native) if native.exists() else None,
        input_manifest_sha256=sha(bundle/'inputs-manifest.json'))
    plan = dict(identity=identity, jobs=all_jobs, status='plan only; no histories', workers=a.workers,
        started_utc=datetime.now(timezone.utc).isoformat(), host=platform.platform(),
        dependencies={name: version(name) for name in ['numpy', 'scipy', 'GitPython']}, python=sys.version,
        histories_per_scenario_program=136000, total_reserved_with_prior=200000)
    if a.resume:
        previous = json.loads((out/'jobs.json').read_text())
        if previous['identity'] != identity or previous['jobs'] != all_jobs:
            raise ValueError('Resume source, inputs, protocol or binary changed')
        plan = previous
    else:
        out.mkdir(parents=True, exist_ok=False)
    def save():
        temp = out/'jobs.json.partial'; temp.write_text(json.dumps(plan, indent=2)+'\n'); temp.replace(out/'jobs.json')
    save()
    if a.plan_only:
        print(f'Frozen {len(all_jobs)} jobs; 136000 new histories/scenario/program; no execution')
        return
    def execute(job):
        target = out/(job['id']+'.tsv'); marker = out/(job['id']+'.done.json')
        if marker.exists():
            done = json.loads(marker.read_text())
            if done['job'] != job or sha(target) != done['output_sha256']:
                raise ValueError('Changed completed batch')
            return job['id']
        argv, env = command(job, protocol, bundle, references, native, probe, target)
        start = time.monotonic()
        recovered = target.exists()
        if not recovered:
            with (out/(job['id']+'.log')).open('w') as log:
                subprocess.run(argv, check=True, stdout=log, stderr=log, env=env)
        # Recover only a fully completed unmarked output. A partial batch is
        # never rerun automatically: its consumed history count is uncertain.
        values = summary(target, job, protocol)
        done = dict(job=job, command=argv, observables=values, output_sha256=sha(target),
                    seconds=time.monotonic()-start, recovered_complete_output=recovered)
        temp = marker.with_suffix('.partial'); temp.write_text(json.dumps(done, indent=2)+'\n'); temp.replace(marker)
        return job['id']
    try:
        for script in ['test_transport_photon_acceptance.py', 'test_transport_photon_pilots.py', 'test_transport_photon_tables.py']:
            with (out/(script+'.log')).open('w') as log:
                subprocess.run([sys.executable, str(ROOT/'tests/crpropa'/script)], check=True,
                    cwd=ROOT, env=ENV, stdout=log, stderr=log)
        plan['status'] = 'checking current production coefficients'; save()
        if not (out/'coefficient-controls.json').exists():
            with ThreadPoolExecutor(max_workers=a.workers) as pool:
                controls = list(pool.map(lambda v: check_coefficients(v, protocol, cache, bundle, references, out), protocol['coefficient_controls']))
            (out/'coefficient-controls.json').write_text(json.dumps(controls, indent=2)+'\n')
        else:
            saved_controls = json.loads((out/'coefficient-controls.json').read_text())
            if [v['control'] for v in saved_controls] != protocol['coefficient_controls']:
                raise ValueError('Coefficient verification inventory changed')
            for control in saved_controls:
                path = out/f'coefficient-{control["control"]["background"]}-{control["control"]["node"]:03}.tsv'
                if sha(path) != control['output_sha256']:
                    raise ValueError('Changed coefficient verification')
        reservation = dict(identity=identity, output=str(out),
            new_histories_per_scenario_program=136000, reserved_total_with_prior=200000)
        ledger = a.campaign_ledger.resolve()
        reserve_campaign(ledger, reservation)
        (out/'campaign-reservation.json').write_text(json.dumps(reservation, indent=2)+'\n')
        plan['status'] = 'running fixed acceptance batches'; save()
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures = [pool.submit(execute, job) for job in all_jobs]
            try:
                for i, future in enumerate(as_completed(futures), 1):
                    print(f'[{i}/{len(all_jobs)}] {future.result()}', flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        plan['status'] = 'execution complete; inspect acceptance-report.json, not wrapper passed'
    except BaseException as error:
        plan.update(status='failed or incomplete; no acceptance', error=str(error))
        raise
    finally:
        save(); report(out, protocol, all_jobs)


if __name__ == '__main__':
    main()
