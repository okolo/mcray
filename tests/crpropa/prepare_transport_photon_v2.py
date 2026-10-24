#!/usr/bin/env python3
"""Prepare fixed ordinary v2 batches from verified previous runtime records.

This performs no nuclear transport and cannot issue acceptance. The complete v2
freeze still requires rare-method pilots, diagnostics, budgets and identities.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_transport_photon_acceptance import PROTOCOL, ROOT, jobs as v1_jobs
from run_transport_photon_pilots import sha

RARE = tuple(prefix+species for species in
             ['electron', 'positron', 'neutrinoE', 'neutrinoAE']
             for prefix in ['yield:', 'energy:'])
V1_SHA = 'e5f345a2773745790f7c59e6941b143cb70450e80a42af38a28aac11deef76a9'
PREVIOUS_COMMIT = '172f25344fff97beed2e41941c32403c46adff53'
PREVIOUS_ARCHIVE_SHA = '75216d7cdaf110500d40c4ad2345410298e03f318e5ef4626c965b0492f0a948'
SEED_BASE = 3000000000


def read(path):
    return json.loads(path.read_text())


def timings(results, sidecar):
    """Verify the completed payload before using execution times only."""
    manifest, outer = read(results/'manifest.json'), read(sidecar)
    if (manifest != {k: v for k, v in outer.items() if k not in ['archive', 'archive_sha256']} or
            manifest['commit'] != PREVIOUS_COMMIT or manifest['status'] != 'passed' or
            manifest['exit_code'] != 0 or outer['archive_sha256'] != PREVIOUS_ARCHIVE_SHA):
        raise ValueError('Wrong or incomplete previous campaign')
    for name, digest in manifest['payload_sha256'].items():
        path = results/name
        if not path.resolve().is_relative_to(results.resolve()) or sha(path) != digest:
            raise ValueError('Unsafe or changed previous payload: '+name)
    protocol = read(PROTOCOL)
    out = results/'transport-photon-acceptance'
    previous = read(out/'jobs.json')
    expected = v1_jobs(protocol)
    if (previous['jobs'] != expected or previous['identity']['protocol_sha256'] != V1_SHA or
            previous['identity']['source_commit'] != PREVIOUS_COMMIT):
        raise ValueError('Previous protocol or job inventory changed')
    samples = defaultdict(list)
    for job in expected:
        done = read(out/(job['id']+'.done.json'))
        if (done['job'] != job or done['recovered_complete_output'] or
                not math.isfinite(done['seconds']) or done['seconds'] <= 0 or
                sha(out/(job['id']+'.tsv')) != done['output_sha256']):
            raise ValueError('Unusable completed timing record: '+job['id'])
        samples[job['case']['background'], job['code'], job['level']].append(
            done['seconds']/job['histories'])
    return {key: dict(seconds_per_primary=statistics.mean(v), batches=len(v),
                     minimum=min(v), maximum=max(v)) for key, v in samples.items()}


def ordinary_plan(protocol, measured, hours=20., workers=16, margin=1.15):
    if (not math.isfinite(hours) or hours <= 0 or isinstance(workers, bool) or
            not isinstance(workers, int) or workers < 1 or
            not math.isfinite(margin) or margin < 1):
        raise ValueError('Invalid runtime budget, workers or timing margin')
    required = {(bg, code, level) for bg in protocol['backgrounds']
                for code in ['mcray', 'crpropa'] for level in range(3)}
    if set(measured) != required:
        raise ValueError('Incomplete timing matrix')
    if any(not math.isfinite(v['seconds_per_primary']) or v['seconds_per_primary'] <= 0
           for v in measured.values()):
        raise ValueError('Invalid timing measurements')
    result, cells = [], []
    for bg_index, bg in enumerate(protocol['backgrounds']):
        for level, (fraction, batches) in enumerate(zip([.5, .25, .25], [128, 64, 64])):
            # Match program histories instead of wasting equal CPU allocations
            # on the much faster official CRPropa adapter. Equal physical
            # workloads also make program comparisons easier to interpret.
            pair_cost = math.fsum(measured[bg, code, level]['seconds_per_primary']
                                  for code in ['mcray', 'crpropa'])
            cpu_budget = hours*3600*workers*fraction/len(protocol['backgrounds'])
            per_batch = math.floor(cpu_budget/(margin*pair_cost*batches))
            if per_batch < 2:
                raise ValueError('Budget too small for independent batches')
            cells.append(dict(background=bg, level=level, batches=batches,
                              histories_per_batch=per_batch,
                              histories_per_program=batches*per_batch,
                              expected_pair_cpu_seconds=batches*per_batch*pair_cost))
            for code_index, code in enumerate(['mcray', 'crpropa']):
                case = dict(id=bg, background=bg, process='talys-decay',
                    pd_photons=True, decay_photons=True,
                    length_mpc=protocol['distance_mpc'],
                    energy_per_nucleon_eV=protocol['energy_per_nucleon_eV'])
                for batch in range(batches):
                    result.append(dict(id=f'v2-{bg}-{code}-l{level}-b{batch:03}',
                        case=case, code=code, level=level, batch=batch,
                        histories=per_batch,
                        seed=SEED_BASE+1000000*code_index+10000*bg_index+batch))
    keys = list(protocol['targets']['cmb/mcray']['means'])
    ordinary = [k for k in keys if k not in RARE]
    if len(ordinary) != 24 or len(set(keys) & set(RARE)) != 8:
        raise ValueError('Changed observable family')
    return dict(jobs=result, cells=cells, ordinary_observables=ordinary,
                rare_observables=list(RARE), total_primary_executions=sum(j['histories'] for j in result),
                expected_ordinary_wall_hours=sum(c['expected_pair_cpu_seconds'] for c in cells)/workers/3600,
                target_wall_hours=hours, workers=workers, timing_margin=margin,
                allocation_by_level=[.5, .25, .25],
                failure_probability_primary=.008, failure_probability_method_checks=.002,
                primary_family_size=384)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--previous-results', required=True, type=Path)
    p.add_argument('--previous-summary', required=True, type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--hours', type=float, default=20.)
    p.add_argument('--workers', type=int, default=16)
    p.add_argument('--timing-margin', type=float, default=1.15)
    a = p.parse_args()
    if a.output is None:
        if 'MCRAY_RESULTS_DIR' not in os.environ:
            p.error('Supply --output or MCRAY_RESULTS_DIR')
        a.output = Path(os.environ['MCRAY_RESULTS_DIR'])/'v2-preparation'
    if sha(PROTOCOL) != V1_SHA:
        raise ValueError('Immutable v1 reference contract changed')
    for test in ['test_transport_photon_v2.py', 'test_rare_channel_sampling.py']:
        subprocess.run([sys.executable, str(Path(__file__).parent/test)], check=True)
    measured = timings(a.previous_results.resolve(), a.previous_summary.resolve())
    plan = ordinary_plan(read(PROTOCOL), measured, a.hours, a.workers, a.timing_margin)
    plan.update(schema='transport-photon-v2-preparation-1',
        status='preparation only; complete v2 freeze and execution adapter pending',
        acceptance_ready=False, created_utc=datetime.now(timezone.utc).isoformat(),
        v1_protocol_sha256=V1_SHA, previous_archive_sha256=PREVIOUS_ARCHIVE_SHA,
        previous_summary_sha256=sha(a.previous_summary),
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        source_diff_sha256=hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD'], cwd=ROOT)).hexdigest(),
        preparer_sha256=sha(Path(__file__)),
        timings={'/'.join(map(str, key)): value for key, value in sorted(measured.items())},
        prior_reserved_per_scenario_program=200000,
        rare_budget=None, rare_proposal=None, method_checks=None,
        note='No new histories. No execution command. Runtime estimate is not a stopping rule. '
             'Fresh seeds; paired numerical settings; no pooling with v1. '
             'Final source/binary/input identities and rare-method freeze are pending.')
    a.output.mkdir(parents=True, exist_ok=False)
    path = a.output/'ordinary-preparation.json'
    path.write_text(json.dumps(plan, indent=2, allow_nan=False)+'\n')
    path.with_suffix('.sha256').write_text(sha(path)+'  '+path.name+'\n')
    print(f'Prepared {len(plan["jobs"])} jobs; {plan["total_primary_executions"]} planned primary executions; no nuclear histories run')


if __name__ == '__main__':
    main()
