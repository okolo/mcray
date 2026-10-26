#!/usr/bin/env python3
"""One fixed, time-limited ordinary transport check before a limited release."""
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

from prepare_transport_photon_v2 import ordinary_plan, timings, RARE, V1_SHA
from run_transport_photon_acceptance import (ROOT, PROTOCOL, ENV, REFERENCE_LIBRARY_SHA256,
    command, summary, report, reserve_campaign, check_coefficients)
from run_transport_photon_pilots import sha, PIN
from run_rare_transport_pilots import begin_batch, verify_inputs

SEED_BASE = 3300000000


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    temporary = path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def make_plan(protocol, measured, workers):
    # Reserve 30 minutes for preparation. Runtime margin leaves additional slack.
    plan = ordinary_plan(protocol, measured, hours=7.5, workers=workers, margin=1.25)
    for job in plan['jobs']:
        job['id'] = job['id'].replace('v2-', 'release8h-', 1)
        job['seed'] += SEED_BASE-3000000000
    # Balance backgrounds/programs/resolutions independently of observed results.
    plan['jobs'].sort(key=lambda j: (j['batch'], j['level'], j['case']['background'], j['code']))
    plan.update(schema='limited-release-check-8h-v1', wall_limit_hours=8,
        release_policy='Release after this check with documented unresolved validation limits; '
                       'no new rare-method campaign and no claim of complete physical acceptance',
        sampling='ordinary production transport only; no importance sampling or splitting',
        main_batches=128, refinement_batches=64)
    return plan


def finish_report(out, protocol, jobs):
    result = report(out, protocol, jobs)
    for row in result['observables']:
        if row['observable'] in RARE:
            row['diagnostic_status'] = row['status']
            row['status'] = 'inconclusive'
            row['role'] = 'rare diagnostic only; interval applicability unestablished'
        else:
            row['role'] = 'ordinary comparison at unchanged tolerance'
    ordinary = [r for r in result['observables'] if r['observable'] not in RARE]
    status = 'inconclusive'
    if any(r['status'] == 'disagreement' for r in ordinary):
        status = 'disagreement'
    elif result['execution_complete'] and all(r['status'] == 'agreement' for r in ordinary):
        status = 'agreement'
    result.update(status='limited release check; not full physical acceptance',
        ordinary_comparison_status=status, full_physical_acceptance=False,
        main_family_status='see ordinary_comparison_status; rare family unresolved',
        confidence_interpretation='Nominal 0.008 Bonferroni allocation over all 384 rows; '
            'rare intervals are exploratory and cannot establish acceptance',
        completed_jobs=sum((out/(j['id']+'.done.json')).exists() for j in jobs),
        expected_jobs=len(jobs), refinement_diagnostics=[])
    if (out/'effective-protocol.json').exists():
        result['effective_protocol_sha256'] = sha(out/'effective-protocol.json')
    for bg in protocol['backgrounds']:
        for code in ['mcray', 'crpropa']:
            by_level = {level: {} for level in range(3)}
            for j in jobs:
                marker = out/(j['id']+'.done.json')
                if j['case']['background'] == bg and j['code'] == code and marker.exists():
                    by_level[j['level']][j['batch']] = read(marker)['observables']
            for level in [1, 2]:
                pairs = sorted(by_level[0].keys() & by_level[level].keys())
                if len(pairs) < 2:
                    continue
                for key in protocol['targets'][bg+'/'+code]['means']:
                    differences = [by_level[level][b][key]-by_level[0][b][key] for b in pairs]
                    result['refinement_diagnostics'].append(dict(background=bg, code=code,
                        level=level, observable=key, paired_batches=len(pairs),
                        mean_difference=statistics.mean(differences),
                        paired_standard_error=statistics.stdev(differences)/math.sqrt(len(pairs)),
                        role='paired numerical diagnostic, not a correction or automatic convergence verdict'))
    result['note'] += ' Release authorization does not convert inconclusive evidence into agreement. '
    write(out/'acceptance-report.json', result)
    write(out/'release-check-report.json', result)
    return result


def execute(job, protocol, bundle, references, native, probe, out, deadline):
    target = out/(job['id']+'.tsv')
    argv, env = command(job, protocol, bundle, references, native, probe, target)
    if time.time() >= deadline:
        return False
    # No automatic replay or recovery in this one-shot campaign.
    if begin_batch(out, job, argv):
        raise ValueError('Started output already exists; this one-shot runner never replays histories')
    start = time.monotonic()
    remaining = deadline-time.time()
    if remaining <= 0:
        raise TimeoutError('Budget expired after start marker; audit before any retry')
    with (out/(job['id']+'.log')).open('w') as log:
        subprocess.run(argv, check=True, stdout=log, stderr=log, env=env, timeout=remaining)
    values = summary(target, job, protocol)
    write(out/(job['id']+'.done.json'), dict(job=job, command=argv, observables=values,
        output_sha256=sha(target), started_sha256=sha(out/(job['id']+'.started.json')),
        seconds=time.monotonic()-start, recovered_complete_output=False))
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['bundle', 'references', 'deterministic-cache', 'previous-results',
                 'previous-summary', 'campaign-ledger']:
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--external-server', action='store_true')
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    if not args.plan_only and (not args.external_server or platform.system() != 'Linux' or platform.machine() != 'x86_64'):
        parser.error('Requires an external native Linux amd64 server')
    if sha(PROTOCOL) != V1_SHA:
        raise ValueError('Frozen physical reference contract changed')
    protocol = read(PROTOCOL)
    measured = timings(args.previous_results.resolve(), args.previous_summary.resolve())
    plan = make_plan(protocol, measured, args.workers)
    out = Path(os.environ['MCRAY_RESULTS_DIR'])/'release-check-8h'
    out.mkdir(parents=True, exist_ok=False)
    plan.update(source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        physical_protocol_sha256=V1_SHA, previous_summary_sha256=sha(args.previous_summary),
        prior_v1_reserved_per_scenario_program=200000,
        python=sys.version, host=platform.platform(),
        dependencies={name: version(name) for name in ['numpy', 'scipy', 'GitPython']},
        created_utc=datetime.now(timezone.utc).isoformat(), status='plan only')
    write(out/'jobs.json', plan)
    if args.plan_only:
        return
    if subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT, text=True).strip():
        raise ValueError('Requires clean frozen source')
    bundle, references, cache = args.bundle.resolve(), args.references.resolve(), args.deterministic_cache.resolve()
    verify_inputs(bundle, references, protocol)
    for name, digest in protocol['reference_files_sha256'].items():
        if sha(cache/name) != digest:
            raise ValueError('Changed deterministic cache: '+name)
    native = ROOT/'build/Release/tests/nuclear_scenario'
    probe = references/'reference-build/crpropa_decay_probe'
    library = references/'reference-build/reference/libcrpropa.so'
    build = read(references/'reference-build.json')
    if (sha(probe) != protocol['crpropa_probe_sha256'] or sha(library) != REFERENCE_LIBRARY_SHA256 or
            build['commit'] != PIN or build['probe_sha256'] != sha(probe)):
        raise ValueError('Retained official runtime/provenance changed')
    # The deadline is established before clone/build/test by the external wrapper.
    deadline = float(os.environ['MCRAY_RELEASE_DEADLINE'])-180
    if not math.isfinite(deadline) or deadline <= time.time():
        raise ValueError('Eight-hour deadline absent or expired')
    plan.update(deadline_unix=deadline,
        binary_sha256={str(p): sha(p) for p in [native, probe, library]},
        inputs_manifest_sha256=sha(bundle/'inputs-manifest.json'),
        reference_manifest_sha256=sha(references/'reference-build.json'),
        driver_sha256=sha(Path(__file__)), status='reserved; checking production coefficients')
    write(out/'jobs.json', plan)
    reserve_campaign(args.campaign_ledger.resolve(), dict(output=str(out), plan=plan))
    write(out/'campaign-reservation.json', dict(output=str(out), plan=plan))
    # The immutable v1 tables and tolerances remain untouched. Only batch counts
    # and the already approved v2 primary confidence allocation change here.
    protocol = dict(protocol, main_batches=128, refinement_batches=64, simultaneous_confidence=.992)
    write(out/'effective-protocol.json', protocol)
    plan['effective_protocol_sha256'] = sha(out/'effective-protocol.json')
    try:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            controls = list(pool.map(lambda e: check_coefficients(e, protocol, cache, bundle, references, out),
                                     protocol['coefficient_controls']))
        write(out/'coefficient-controls.json', controls)
        plan['status'] = 'running fixed ordinary batches'
        write(out/'jobs.json', plan)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(execute, j, protocol, bundle, references, native, probe, out, deadline)
                       for j in plan['jobs']]
            try:
                for i, future in enumerate(as_completed(futures), 1):
                    if future.result():
                        print(f'[{i}/{len(futures)}] completed', flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        plan['status'] = 'complete' if all((out/(j['id']+'.done.json')).exists() for j in plan['jobs']) else 'budget exhausted; incomplete'
    except BaseException as error:
        plan.update(status='failed or budget exhausted; inspect retained partial outputs', error=str(error))
        raise
    finally:
        write(out/'jobs.json', plan)
        finish_report(out, protocol, plan['jobs'])
    if plan['status'] != 'complete':
        raise TimeoutError('Fixed campaign incomplete within time budget; do not rerun')


if __name__ == '__main__':
    main()
