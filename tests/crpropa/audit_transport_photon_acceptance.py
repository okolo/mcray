#!/usr/bin/env python3
"""Verify a frozen acceptance archive and recompute decisions from raw histories.

This reads completed outputs only. It generates no histories, changes no
acceptance criteria and never overwrites the archived report.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

from run_transport_photon_acceptance import PROTOCOL, ROOT, jobs, summary, sha


def read(path):
    return json.loads(path.read_text())


def equal_number(actual, expected, scale=0.):
    # Verification guard for Python/SciPy accumulation/quantile roundoff only;
    # this is not an additional physical acceptance tolerance.
    if not math.isfinite(actual) or not math.isfinite(expected):
        raise ValueError('Nonfinite audit value')
    if abs(actual-expected) > 5e-12*max(abs(actual), abs(expected), abs(scale)):
        raise ValueError(f'Numerical report mismatch: {actual}, {expected}')


def audit(results, sidecar, archive, output):
    from scipy.stats import t
    meta, outer = read(results/'manifest.json'), read(sidecar)
    if meta != {k: v for k, v in outer.items() if k not in ['archive', 'archive_sha256']}:
        raise ValueError('Archive manifest and sidecar differ')
    if sha(archive) != outer['archive_sha256']:
        raise ValueError('Archive checksum mismatch')
    if meta['status'] != 'passed' or meta['exit_code'] != 0:
        raise ValueError('Execution did not complete successfully')
    commit = meta['commit']
    if not re.fullmatch('[0-9a-f]{40}', commit):
        raise ValueError('Invalid frozen source revision')
    source_files = ['run_transport_photon_acceptance.py',
                    'run_transport_photon_pilots.py', 'run_network_rate_scan.py']
    for name in source_files+['protocols/transport_photons_v1.json']:
        path = 'tests/crpropa/'+name
        frozen = subprocess.check_output(['git', 'show', commit+':'+path], cwd=ROOT)
        if frozen != (ROOT/path).read_bytes():
            raise ValueError('Use audit helpers matching the frozen source: '+path)
    for name, digest in meta['payload_sha256'].items():
        path = results/name
        if not path.resolve().is_relative_to(results) or sha(path) != digest:
            raise ValueError('Changed or unsafe payload: '+name)
    expected_files = set(meta['payload_sha256']) | {'manifest.json'}
    actual_files = {str(p.relative_to(results)) for p in results.rglob('*') if p.is_file()}
    if actual_files != expected_files:
        raise ValueError('Unexpected result inventory')

    out = results/'transport-photon-acceptance'
    plan, protocol = read(out/'jobs.json'), read(PROTOCOL)
    frozen_jobs = jobs(protocol)
    if plan['jobs'] != frozen_jobs or plan['identity']['source_commit'] != commit:
        raise ValueError('Frozen jobs or source changed')
    identity = plan['identity']
    if (identity['protocol_sha256'] != sha(PROTOCOL) or
            identity['driver_sha256'] != sha(ROOT/'tests/crpropa/run_transport_photon_acceptance.py') or
            identity['source_patch_sha256'] != hashlib.sha256(b'').hexdigest()):
        raise ValueError('Protocol, driver or clean-source identity mismatch')
    reservation = read(out/'campaign-reservation.json')
    if reservation['identity'] != identity or reservation['reserved_total_with_prior'] != protocol['history_cap']:
        raise ValueError('Campaign reservation mismatch')
    controls = read(out/'coefficient-controls.json')
    if [c['control'] for c in controls] != protocol['coefficient_controls']:
        raise ValueError('Coefficient-control inventory mismatch')
    byte_identical = 0
    for c in controls:
        entry = c['control']
        path = out/f'coefficient-{entry["background"]}-{entry["node"]:03}.tsv'
        if c['status'] != 'matched' or sha(path) != c['output_sha256']:
            raise ValueError('Changed coefficient control')
        byte_identical += sha(path) == entry['native_sha256']

    grouped, events, executions = defaultdict(list), defaultdict(Counter), Counter()
    for job in frozen_jobs:
        marker = read(out/(job['id']+'.done.json'))
        path = out/(job['id']+'.tsv')
        if marker['job'] != job or sha(path) != marker['output_sha256']:
            raise ValueError('Changed completed batch: '+job['id'])
        actual = summary(path, job, protocol)
        if actual.keys() != marker['observables'].keys():
            raise ValueError('Changed observable inventory')
        for key, value in actual.items():
            equal_number(value, marker['observables'][key])
        group = (job['case']['background'], job['code'], job['level'])
        grouped[group].append(actual)
        executions[group[:2]] += job['histories']
        for key in ['yield:electron', 'yield:positron', 'yield:photon']:
            events[group][key] += round(actual[key]*job['histories'])
    if any(n != 136000 for n in executions.values()) or len(executions) != 4:
        raise ValueError('History accounting mismatch')

    archived = read(out/'acceptance-report.json')
    family = 3*2*len(protocol['backgrounds'])*len(protocol['targets']['cmb/mcray']['means'])
    if family != archived['simultaneous_family_size'] or family != len(archived['observables']):
        raise ValueError('Incomplete comparison family')
    rows = {(r['background'], r['code'], r['level'], r['observable']): r for r in archived['observables']}
    counts, group_counts, verified = Counter(), {}, []
    for group, batches in sorted(grouped.items()):
        bg, code, level = group
        n = protocol['main_batches'] if level == 0 else protocol['refinement_batches']
        if len(batches) != n:
            raise ValueError('Incomplete batch group')
        target = protocol['targets'][bg+'/'+code]
        normal = protocol['targets'][bg+'/independent']['means']
        critical = float(t.ppf(1-(1-protocol['simultaneous_confidence'])/(2*family), n-1))
        statuses = Counter()
        for key, reference in target['means'].items():
            values = [b[key] for b in batches]
            mean = math.fsum(values)/n
            se = math.sqrt(math.fsum((v-mean)**2 for v in values)/(n-1)/n)
            half, norm = critical*se, abs(normal[key])
            eps = target['numerical_envelope'][key]
            spectral = key.startswith('photon_bin_')
            delta = norm*(protocol['spectral_tolerance'] if spectral else protocol['integral_tolerance'])
            difference, status = mean-reference, 'inconclusive'
            if se > 0 and norm > 0 and math.isfinite(half):
                if abs(difference)-half-eps > delta:
                    status = 'disagreement'
                elif abs(difference)+half+eps <= delta and (spectral or se <= .01*norm):
                    status = 'agreement'
            original = rows[(*group, key)]
            values_to_check = dict(mean=mean, target=reference, standard_error=se,
                                   half_width=half, normalization=norm, tolerance=delta,
                                   numerical_envelope=eps, difference=difference)
            for field, value in values_to_check.items():
                equal_number(value, original[field], norm)
            if status != original['status']:
                raise ValueError('Decision changed: '+str((*group, key)))
            statuses[status] += 1
            verified.append(dict(original))
        counts.update(statuses)
        group_counts['/'.join(map(str, group))] = dict(statuses)
    coefficient_rows = archived['deterministic_coefficient_comparisons']
    for row in coefficient_rows:
        bg, code, key = row['background'], row['code'], row['observable']
        ref = protocol['targets'][bg+'/independent']
        native = protocol['targets'][bg+'/'+code]
        norm = abs(ref['means'][key])
        eps = ref['numerical_envelope'][key]+native['numerical_envelope'][key]
        diff = native['means'][key]-ref['means'][key]
        delta = norm*(protocol['spectral_tolerance'] if key.startswith('photon_bin_') else protocol['integral_tolerance'])
        status = ('agreement' if norm and abs(diff)+eps <= delta else
                  'disagreement' if norm and abs(diff)-eps > delta else 'inconclusive')
        if row['status'] != status:
            raise ValueError('Coefficient decision changed')
    native_rows = [r for r in coefficient_rows if r['code'] == 'mcray']
    overall = ('disagreement' if counts['disagreement'] or any(r['status'] == 'disagreement' for r in native_rows) else
               'agreement' if counts['agreement'] == family and all(r['status'] == 'agreement' for r in native_rows)
               else 'inconclusive')
    if archived['status'] != overall or not archived['execution_complete']:
        raise ValueError('Global verdict changed or incomplete')
    ctest = ET.parse(results/'ctest.xml').getroot()
    if ctest.get('failures') != '0' or ctest.get('skipped') != '0':
        raise ValueError('Prerequisite checks failed or were skipped')
    result = dict(status='audited completed frozen campaign', physical_status=overall,
        source_commit=commit, archive_sha256=outer['archive_sha256'],
        protocol_sha256=sha(PROTOCOL), payload_files_verified=len(meta['payload_sha256']),
        jobs_verified=len(frozen_jobs), new_primary_executions=sum(executions.values()),
        histories_per_scenario_program={'/'.join(k): v for k, v in executions.items()},
        comparisons=dict(counts), groups=group_counts, coefficient_controls=len(controls),
        byte_identical_coefficient_controls=byte_identical, ctest=ctest.attrib,
        native_coefficient_statuses=dict(Counter(r['status'] for r in native_rows)),
        observed_event_counts={'/'.join(map(str, k)): dict(v) for k, v in events.items()},
        report_rows=verified, dependencies=plan['dependencies'],
        analysis_python=sys.version, auditor_sha256=sha(Path(__file__)),
        note='No new histories; archived files unchanged; roundoff guard is not a physical tolerance.')
    output.mkdir(parents=True, exist_ok=False)
    (output/'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(overall, dict(counts), 'verified raw batches:', len(frozen_jobs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True, type=Path)
    parser.add_argument('--summary', required=True, type=Path)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    audit(args.results.resolve(), args.summary.resolve(), args.archive.resolve(), args.output.resolve())
