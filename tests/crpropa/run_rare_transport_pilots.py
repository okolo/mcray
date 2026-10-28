#!/usr/bin/env python3
"""Independent full-Fe56 weighted-method pilots. Never issues acceptance."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import version
import json
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from run_transport_photon_acceptance import (PROTOCOL, ROOT, ENV, REFERENCE_LIBRARY_SHA256,
    command as ordinary_command, reserve_campaign)
from run_transport_photon_pilots import PIN, sha
from run_network_rate_scan import verify
from prepare_rare_guides import prepare
from prepare_transport_photon_v2 import V1_SHA
from rare_pilot_statistics import summarize, batch_statistics

SAMPLING_BASE = '54b24e0fa7e6aa55c89804b51afd7611f632945d'


def read(path):
    return json.loads(path.read_text())


def begin_batch(out, job, argv):
    """An absent output does not prove that a crashed native job used no histories."""
    target, started = out/(job['id']+'.tsv'), out/(job['id']+'.started.json')
    if started.exists():
        previous = read(started)
        if previous['job'] != job or previous['command'] != argv:
            raise ValueError('Changed started batch')
        if not target.exists():
            raise ValueError('Started batch has no output; consumed histories unknown, never automatically rerun')
        return True
    if target.exists():
        raise ValueError('Output has no launcher start provenance; audit before recovery')
    with started.open('x') as stream:
        json.dump(dict(job=job, command=argv, started_utc=datetime.now(timezone.utc).isoformat()), stream, indent=2)
        stream.write('\n')
    return False


def jobs(protocol, profile):
    if profile not in ['smoke', 'pilot']:
        raise ValueError('Unknown pilot profile')
    settings = [('ordinary', 'minus', 0., 1), ('split', 'minus', 0., 4),
                ('minus-moderate', 'minus', .5, 4), ('minus-strong', 'minus', .8, 4),
                ('plus-moderate', 'plus', .5, 4), ('plus-strong', 'plus', .8, 4)]
    batches, size, levels = (1, 2, 1) if profile == 'smoke' else (8, 500, 3)
    result = []
    for bi, bg in enumerate(protocol['backgrounds']):
        case = dict(id=bg, background=bg, process='talys-decay', pd_photons=True, decay_photons=True,
                    length_mpc=protocol['distance_mpc'], energy_per_nucleon_eV=protocol['energy_per_nucleon_eV'])
        for ci, code in enumerate(['mcray', 'crpropa']):
            for pi, (name, target, fraction, copies) in enumerate(settings):
                for level in range(levels):
                    for batch in range(batches):
                        result.append(dict(id=f'{bg}-{code}-{name}-l{level}-b{batch:02}', case=case,
                            code=code, profile=name, target=target, fraction=fraction, copies=copies,
                            level=level, batch=batch, histories=size,
                            # New pilot namespace, distinct from v1 and v2 ordinary seeds.
                            # Levels alone are paired; methods/programs/backgrounds are independent.
                            seed=3200000000+(10000000 if profile == 'smoke' else 0)+bi*1000000+ci*500000+pi*50000+batch))
    return result


def verify_inputs(bundle, references, protocol):
    inputs = read(bundle/'inputs-manifest.json')
    if inputs['data_commit'] != '34b10e27375c0e9894d7c984ac15fb96e12975ef' or inputs['reference_code'] != PIN:
        raise ValueError('Changed frozen input versions')
    for name, digest in inputs['payload_sha256'].items():
        path = bundle/name
        if not path.resolve().is_relative_to(bundle) or sha(path) != digest:
            raise ValueError('Changed or unsafe input: '+name)
    for name, digest in protocol['data_manifest_sha256'].items():
        if sha(bundle/'nuclei'/name) != digest:
            raise ValueError('Changed frozen nuclear manifest')
    for bg, digest in protocol['prepared_reference_manifest_sha256'].items():
        verify(references/bg, 'shared_manifest.json')
        if sha(references/bg/'shared_manifest.json') != digest:
            raise ValueError('Changed prepared reference')


def build_reference(out, logdir):
    source, build = out/'reference-source', out/'reference-build'
    marker = out/'reference-build.json'
    if marker.exists():
        meta = read(marker)
        if meta['commit'] != PIN:
            raise ValueError('Changed reference source identity')
        for name, digest in meta['binary_sha256'].items():
            if sha(build/name) != digest:
                raise ValueError('Changed weighted reference binary')
        return build
    if source.exists() or build.exists():
        raise ValueError('Incomplete reference build: retain artifacts; do not reuse a partial campaign')
    commands = [
        ['git', 'clone', '--no-checkout', 'https://github.com/CRPropa/CRPropa3.git', str(source)],
        ['git', '-C', str(source), 'fetch', 'origin', PIN],
        ['git', '-C', str(source), 'checkout', '--detach', 'FETCH_HEAD'],
        ['cmake', '-S', str(ROOT/'tests/crpropa'), '-B', str(build), '-DCRPROPA_SOURCE_DIR='+str(source),
         '-DBUILD_DOCKER_MCRAY_PROBE=OFF', '-DCMAKE_BUILD_TYPE=Release'],
        ['cmake', '--build', str(build), '--target', 'crpropa_weighted_probe', 'crpropa_decay_probe', '--parallel', '2']]
    for i, argv in enumerate(commands):
        with (logdir/f'reference-build-{i}.log').open('w') as log:
            subprocess.run(argv, check=True, env=ENV, stdout=log, stderr=subprocess.STDOUT)
    actual = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if actual != PIN:
        raise ValueError('Wrong weighted reference source revision')
    paths = [build/'crpropa_weighted_probe', build/'crpropa_decay_probe', build/'reference/libcrpropa.so']
    marker.write_text(json.dumps(dict(commit=actual,
        binary_sha256={str(p.relative_to(build)): sha(p) for p in paths}), indent=2)+'\n')
    return build


def command(job, protocol, bundle, references, native_build, ref_build, out, local):
    target = out/(job['id']+'.tsv')
    if job['profile'] == 'ordinary':
        probe = (ref_build if local else references/'reference-build')/'crpropa_decay_probe'
        argv, env = ordinary_command(job, protocol, bundle, references,
            native_build/'nuclear_scenario', probe, target)
        if local and job['code'] == 'crpropa':
            env = {**env, 'DYLD_LIBRARY_PATH': str(ref_build/'reference'), 'LD_LIBRARY_PATH': str(ref_build/'reference')}
        return argv, env
    common = ['--guide', str(out/'guides'/(job['target']+'.tsv')), '--output', str(target),
              '--primaries', str(job['histories']), '--seed', str(job['seed']),
              '--fraction', str(job['fraction']), '--copies', str(job['copies']),
              '--distance-mpc', str(protocol['distance_mpc'])]
    if job['code'] == 'mcray':
        control = protocol['native_controls'][job['level']]
        return ([str(native_build/'weighted_nuclear_scenario'), '--data', str(bundle/'nuclei'),
            '--photon-table', str(references/job['case']['background']/'shared_photons.tsv'),
            '--energy-per-nucleon', str(protocol['energy_per_nucleon_eV']),
            '--background-step', str(control['step']), '--eps', str(control['eps'])]+common, ENV)
    masses = {tuple(map(int, r[:2])): float(r[2]) for line in (bundle/'nuclei/masses.tsv').read_text().splitlines()
              if line and not line.startswith('#') for r in [line.split()]}
    gamma = 56*protocol['energy_per_nucleon_eV']/(masses[56, 26]*1e6)
    env = {**ENV, 'CRPROPA_DATA_PATH': str(references/job['case']['background']),
           'LD_LIBRARY_PATH': str(ref_build/'reference'), 'DYLD_LIBRARY_PATH': str(ref_build/'reference')}
    return ([str(ref_build/'crpropa_weighted_probe'), '--gamma', str(gamma),
             '--step-fraction', str(protocol['reference_step_fractions'][job['level']])]+common, env)


def report(out, all_jobs):
    groups, complete = defaultdict(list), 0
    total_extra = 0
    for job in all_jobs:
        marker = out/(job['id']+'.done.json')
        if not marker.exists():
            continue
        done = read(marker)
        if done['job'] != job or sha(out/(job['id']+'.tsv')) != done['output_sha256']:
            raise ValueError('Changed completed pilot batch')
        complete += 1
        total_extra += done['summary']['extra_continuations']
        groups[job['case']['background'], job['code'], job['level'], job['profile']].append(done)
    summaries = {'/'.join(map(str, key)): batch_statistics(value) for key, value in sorted(groups.items())}
    efficiency = []
    for key, value in groups.items():
        bg, code, level, profile = key
        baseline = summaries.get(f'{bg}/{code}/{level}/ordinary', {})
        stats = summaries['/'.join(map(str, key))]
        if profile == 'ordinary' or 'observables' not in baseline or 'observables' not in stats:
            continue
        for observable, estimate in stats['observables'].items():
            v, v0 = estimate['batch_variance'], baseline['observables'][observable]['batch_variance']
            efficiency.append(dict(background=bg, code=code, level=level, profile=profile, observable=observable,
                variance_ratio=v0/v if v > 0 and v0 > 0 else None,
                time_adjusted_gain=(v0*baseline['mean_batch_seconds']/(v*stats['mean_batch_seconds'])) if
                    v > 0 and v0 > 0 and baseline['mean_batch_seconds'] is not None and stats['mean_batch_seconds'] is not None else None,
                note='Unresolved zero variance gives no efficiency claim; pilot only'))
    result = dict(status='independent method pilots only; no physical acceptance',
        completed_jobs=complete, expected_jobs=len(all_jobs), extra_continuations=total_extra,
        completed_primary_executions=sum(j['histories'] for j in all_jobs if (out/(j['id']+'.done.json')).exists()),
        unmarked_outputs=[j['id'] for j in all_jobs if (out/(j['id']+'.tsv')).exists() and
                          not (out/(j['id']+'.done.json')).exists()],
        incomplete_started_batches=[j['id'] for j in all_jobs if (out/(j['id']+'.started.json')).exists() and
                                    not (out/(j['id']+'.done.json')).exists()],
        preflight_status='complete' if all((out/name/'controls.json').exists() for name in ['identity-controls', 'decay-controls'])
            else 'incomplete or not started; inspect outputs before any retry',
        groups=summaries, efficiency=efficiency,
        pending=['pilot analysis of weight tails and batch interval applicability',
                 'frozen rare proposal, sample budget and statistical validation inventory',
                 'new final acceptance seeds and identities; no pooling these pilots'])
    (out/'pilot-report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle', required=True, type=Path)
    p.add_argument('--references', required=True, type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--native-build', type=Path, default=ROOT/'build/Release/tests')
    p.add_argument('--reference-build', type=Path, help='Development smoke only: prebuilt pinned reference')
    p.add_argument('--campaign-ledger', type=Path)
    p.add_argument('--workers', type=int, default=16)
    p.add_argument('--profile', choices=['smoke', 'pilot'], default='pilot')
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--external-server', action='store_true')
    p.add_argument('--local-smoke', action='store_true')
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    if a.workers < 1:
        p.error('Workers must be positive')
    if a.local_smoke and (a.profile != 'smoke' or a.reference_build is None):
        p.error('Local execution is limited to smoke with an explicit prebuilt reference')
    if not a.plan_only and not a.local_smoke and (not a.external_server or platform.system() != 'Linux' or platform.machine() != 'x86_64'):
        p.error('Long pilots require an external native Linux amd64 server')
    if a.reference_build is not None and not a.local_smoke:
        p.error('Production pilots build their own pinned reference; preserve the retained official control')
    if not a.plan_only and a.campaign_ledger is None:
        p.error('Execution requires a persistent campaign ledger')
    if a.output is None:
        if 'MCRAY_RESULTS_DIR' not in os.environ:
            p.error('Supply --output or MCRAY_RESULTS_DIR')
        a.output = Path(os.environ['MCRAY_RESULTS_DIR'])/'rare-transport-pilots'
    bundle, references, out = a.bundle.resolve(), a.references.resolve(), a.output.resolve()
    if sha(PROTOCOL) != V1_SHA:
        raise ValueError('Immutable v1 physics/reference contract changed')
    protocol = read(PROTOCOL)
    verify_inputs(bundle, references, protocol)
    subprocess.run(['git', 'merge-base', '--is-ancestor', SAMPLING_BASE, 'HEAD'], cwd=ROOT, check=True)
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT, text=True).strip()
    if dirty and not a.local_smoke and not a.plan_only:
        p.error('Server execution requires a clean published source commit')
    all_jobs = jobs(protocol, a.profile)
    source_files = ['run_rare_transport_pilots.py', 'prepare_rare_guides.py', 'rare_pilot_statistics.py',
                    'check_weighted_adapters.py', 'check_weighted_decay.py',
                    'weighted_probe.cpp', 'weighted_common.cpp', 'weighted_common.h']
    identity = dict(source_commit=source, protocol_sha256=sha(PROTOCOL), sampling_fix_ancestor=SAMPLING_BASE,
        source_patch_sha256=hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD'], cwd=ROOT)).hexdigest(),
        native_adapter_sha256=sha(ROOT/'tests/weighted_nuclear_scenario.cpp'),
        helper_sha256={name: sha(ROOT/'tests/crpropa'/name) for name in source_files},
        bundle_manifest_sha256=sha(bundle/'inputs-manifest.json'),
        reference_manifests={bg: sha(references/bg/'shared_manifest.json') for bg in protocol['backgrounds']},
        profile=a.profile, python=sys.version,
        dependencies={name: version(name) for name in ['numpy', 'scipy', 'GitPython']})
    if a.resume:
        previous = read(out/'jobs.json')
        if previous['identity'] != identity or previous['jobs'] != all_jobs or previous['status'] == 'plan only':
            raise ValueError('Resume identity or fixed jobs changed')
        for name, digest in previous['guide_manifest']['payload_sha256'].items():
            if sha(out/'guides'/name) != digest:
                raise ValueError('Changed importance guide')
        plan = previous
    else:
        out.mkdir(parents=True, exist_ok=False)
        plan = dict(schema='rare-transport-method-pilots-v1', identity=identity, jobs=all_jobs, status='preparing independent pilots',
            created_utc=datetime.now(timezone.utc).isoformat(), workers=a.workers,
            guide_manifest=prepare(bundle/'nuclei', out/'guides'),
            prior_reserved_per_scenario_program=200000,
            new_primary_executions=sum(j['histories'] for j in all_jobs),
            preflight_control_primary_executions=48,
            finite_decay_control_primary_executions=1024,
            max_extra_continuations=sum(j['histories']*(j['copies']-1) for j in all_jobs),
            note='No final acceptance. Additional pilot namespace; no reuse of v1 or v2 final samples.')
        plan['proposal'] = 'defensive mixture of physical probabilities and topology density on positive physical support'
    def save():
        tmp = out/'jobs.json.partial'; tmp.write_text(json.dumps(plan, indent=2)+'\n'); tmp.replace(out/'jobs.json')
    save()
    if a.plan_only:
        plan['status'] = 'plan only'; save(); print(f'{len(all_jobs)} fixed jobs; {plan["new_primary_executions"]} planned primary executions; no histories run'); return
    for test in ['test_rare_transport_pilots.py', 'test_rare_channel_sampling.py']:
        with (out/(test+'.log')).open('w') as log:
            subprocess.run([sys.executable, str(ROOT/'tests/crpropa'/test)], check=True,
                           env=ENV, stdout=log, stderr=subprocess.STDOUT)
    # Do not overwrite or relink the retained official runtime.
    if not a.local_smoke:
        if (sha(references/'reference-build/crpropa_decay_probe') != protocol['crpropa_probe_sha256'] or
                sha(references/'reference-build/reference/libcrpropa.so') != REFERENCE_LIBRARY_SHA256):
            raise ValueError('Changed retained official reference binary/library')
    native_build = a.native_build.resolve()
    ref_build = a.reference_build.resolve() if a.local_smoke else build_reference(out, out)
    binaries = [native_build/'weighted_nuclear_scenario', native_build/'nuclear_scenario',
                ref_build/'crpropa_weighted_probe', ref_build/'crpropa_decay_probe']
    binaries += list((ref_build/'reference').glob('libcrpropa.so'))+list((ref_build/'reference').glob('libcrpropa.dylib'))
    fingerprints = {str(path): sha(path) for path in binaries}
    if a.resume and plan['binary_sha256'] != fingerprints:
        raise ValueError('Resume binaries changed')
    plan['binary_sha256'] = fingerprints
    reservation = dict(identity=identity, output=str(out), jobs=all_jobs,
                       binary_sha256=fingerprints,
                       primary_reservation=plan['new_primary_executions']+plan['preflight_control_primary_executions']+
                           plan['finite_decay_control_primary_executions'])
    reserve_campaign(a.campaign_ledger.resolve(), reservation)
    (out/'campaign-reservation.json').write_text(json.dumps(reservation, indent=2)+'\n')
    plan['status'] = 'running independent pilots'; save()
    def identity_controls():
        controls_path = out/'identity-controls/controls.json'
        if controls_path.exists():
            controls = read(controls_path)
            if controls['primary_executions'] != 48 or len(controls['controls']) != 12:
                raise ValueError('Changed preflight identity inventory')
            for entry in controls['controls']:
                name = entry['job']['id']
                if (sha(out/'identity-controls'/(name+'.tsv')) != entry['output_sha256'] or
                        sha(out/'identity-controls'/(name+'-ordinary.tsv')) != entry['ordinary_sha256']):
                    raise ValueError('Changed preflight output')
        else:
            from check_weighted_adapters import check
            check(bundle, references, native_build, ref_build, out/'guides', out/'identity-controls')
        decay_path = out/'decay-controls/controls.json'
        if decay_path.exists():
            controls = read(decay_path)
            if controls['primary_executions'] != 1024 or len(controls['controls']) != 4:
                raise ValueError('Changed finite-decay control inventory')
            for entry in controls['controls']:
                if sha(out/'decay-controls'/(entry['id']+'.tsv')) != entry['output_sha256']:
                    raise ValueError('Changed finite-decay control output')
        else:
            from check_weighted_decay import check
            check(bundle, references, native_build, ref_build, out/'guides', out/'decay-controls')
    def execute(job):
        target, marker = out/(job['id']+'.tsv'), out/(job['id']+'.done.json')
        if marker.exists():
            done = read(marker)
            if (done['job'] != job or sha(target) != done['output_sha256'] or
                    sha(out/(job['id']+'.started.json')) != done['started_sha256'] or
                    summarize(target, job, protocol) != done['summary']):
                raise ValueError('Changed completed pilot')
            return job['id']
        argv, env = command(job, protocol, bundle, references, native_build, ref_build, out, a.local_smoke)
        recovered = begin_batch(out, job, argv)
        start = time.monotonic()
        if not recovered:
            with (out/(job['id']+'.log')).open('w') as log:
                subprocess.run(argv, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
        stats = summarize(target, job, protocol)
        done = dict(job=job, command=argv, seconds=time.monotonic()-start,
                    started_sha256=sha(out/(job['id']+'.started.json')),
                    recovered_complete_output=recovered, output_sha256=sha(target), summary=stats)
        tmp = marker.with_suffix('.partial'); tmp.write_text(json.dumps(done, indent=2, allow_nan=False)+'\n'); tmp.replace(marker)
        return job['id']
    try:
        plan['status'] = 'checking zero-bias identities'; save()
        identity_controls()
        plan['status'] = 'running independent pilots'; save()
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures = [pool.submit(execute, job) for job in all_jobs]
            try:
                for i, future in enumerate(as_completed(futures), 1):
                    print(f'[{i}/{len(all_jobs)}] {future.result()}', flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        plan['status'] = 'complete independent method pilots; final acceptance freeze pending'
    except BaseException as error:
        plan.update(status='failed or incomplete; no acceptance', error=str(error)); raise
    finally:
        save(); report(out, all_jobs)


if __name__ == '__main__':
    main()
