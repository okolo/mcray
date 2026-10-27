#!/usr/bin/env python3
"""Recompute retained ordinary release-check results without running transport."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import statistics

from prepare_transport_photon_v2 import RARE, V1_SHA
from run_transport_photon_acceptance import PROTOCOL, summary, interval
from run_transport_photon_pilots import sha


def read(path):
    return json.loads(path.read_text())


def close(a, b):
    if a.keys() != b.keys() or any(not math.isclose(a[k], b[k], rel_tol=5e-13, abs_tol=1e-25) for k in a):
        raise ValueError('Stored batch statistics differ from raw histories')


def audit(results, sidecar, output, archive=None):
    outer, manifest = read(sidecar), read(results/'manifest.json')
    if manifest != {k: v for k, v in outer.items() if k not in ['archive', 'archive_sha256']}:
        raise ValueError('Outer and embedded manifests differ')
    if archive is not None and sha(archive) != outer['archive_sha256']:
        raise ValueError('Archive checksum mismatch')
    for name, digest in manifest['payload_sha256'].items():
        path = results/name
        if not path.resolve().is_relative_to(results.resolve()) or sha(path) != digest:
            raise ValueError('Unsafe or changed payload: '+name)
    data = results/'release-check-8h'
    plan, protocol = read(data/'jobs.json'), read(data/'effective-protocol.json')
    if (sha(PROTOCOL) != V1_SHA or plan['physical_protocol_sha256'] != V1_SHA or
            protocol != dict(read(PROTOCOL), main_batches=128, refinement_batches=64, simultaneous_confidence=.992) or
            sha(data/'effective-protocol.json') != plan['effective_protocol_sha256'] or
            plan['source_commit'] != manifest['commit'] or
            plan['schema'] != 'limited-release-check-8h-v1'):
        raise ValueError('Source or frozen reference contract changed')
    controls = read(data/'coefficient-controls.json')
    if [c['control'] for c in controls] != protocol['coefficient_controls']:
        raise ValueError('Coefficient-control inventory changed')
    for c in controls:
        name = f'coefficient-{c["control"]["background"]}-{c["control"]["node"]:03}.tsv'
        if c['status'] != 'matched' or sha(data/name) != c['output_sha256']:
            raise ValueError('Coefficient-control output changed')
    groups, expected, incomplete = defaultdict(dict), Counter(), []
    completed_histories = 0
    for i, job in enumerate(plan['jobs']):
        group = job['case']['background'], job['code'], job['level']
        expected[group] += 1
        marker = data/(job['id']+'.done.json')
        if not marker.exists():
            # A killed native process can have consumed histories before any
            # products reach its output. Never infer zero consumption here.
            incomplete.append(dict(job=job, consumption='unknown; at most reserved histories',
                output_bytes=(data/(job['id']+'.tsv')).stat().st_size if (data/(job['id']+'.tsv')).exists() else None))
            continue
        done = read(marker)
        started = data/(job['id']+'.started.json')
        if (done['job'] != job or sha(started) != done['started_sha256'] or
                read(started)['command'] != done['command'] or
                read(started)['job'] != job or done['recovered_complete_output'] or
                sha(data/(job['id']+'.tsv')) != done['output_sha256']):
            raise ValueError('Completed job provenance changed')
        values = summary(data/(job['id']+'.tsv'), job, protocol)
        close(values, done['observables'])
        groups[group][job['batch']] = dict(values=values, histories=job['histories'])
        completed_histories += job['histories']
        if i % 100 == 0:
            print(f'Audited {i+1}/{len(plan["jobs"])} planned jobs', flush=True)
    rows, refinements, inventory = [], [], []
    for group, count in sorted(expected.items()):
        bg, code, level = group
        batches = groups[group]
        complete = len(batches) == count
        inventory.append(dict(background=bg, code=code, level=level, completed_batches=len(batches),
            planned_batches=count, completed_histories=sum(b['histories'] for b in batches.values())))
        if len(batches) < 2:
            continue
        target = protocol['targets'][bg+'/'+code]
        norms = protocol['targets'][bg+'/independent']['means']
        for key in target['means']:
            estimates = [batches[b]['values'][key] for b in sorted(batches)]
            row = interval(estimates, target['means'][key], abs(norms[key]),
                target['numerical_envelope'][key], key.startswith('photon_bin_'), protocol, 384)
            row.update(background=bg, code=code, level=level, observable=key, batches=len(batches),
                frozen_group_complete=complete, nominal_status=row['status'],
                role='ordinary comparison' if key not in RARE else 'rare diagnostic only')
            if not complete or key in RARE:
                row['status'] = 'inconclusive'
            if key.startswith('yield:'):
                row['observed_particles'] = round(sum(v['values'][key]*v['histories'] for v in batches.values()))
            rows.append(row)
    for (bg, code, level), batches in sorted(groups.items()):
        if level == 0:
            continue
        baseline = groups[bg, code, 0]
        pairs = sorted(batches.keys() & baseline.keys())
        if len(pairs) < 2:
            continue
        for key in protocol['targets'][bg+'/'+code]['means']:
            differences = [batches[b]['values'][key]-baseline[b]['values'][key] for b in pairs]
            norm = abs(protocol['targets'][bg+'/independent']['means'][key])
            refinements.append(dict(background=bg, code=code, level=level, observable=key,
                paired_batches=len(pairs), mean_difference=statistics.mean(differences),
                paired_standard_error=statistics.stdev(differences)/math.sqrt(len(pairs)),
                normalization=norm, role='paired diagnostic; no empirical mean correction'))
    result = dict(schema='retained-release-check-audit-v1', archive=outer['archive'],
        archive_sha256=outer['archive_sha256'], source_commit=manifest['commit'],
        payloads_verified=len(manifest['payload_sha256']), source_status=manifest['status'],
        completed_jobs=sum(len(v) for v in groups.values()), expected_jobs=len(plan['jobs']),
        completed_primary_histories=completed_histories, reserved_primary_histories=plan['total_primary_executions'],
        incomplete_jobs=incomplete, group_inventory=inventory, observables=rows, refinements=refinements,
        status='incomplete resource-limited release check' if incomplete else 'completed resource-limited release check',
        full_physical_acceptance=False,
        note='Offline recomputation only; no histories rerun. Incomplete groups and rare outputs cannot '
             'establish acceptance. Nominal intervals for incomplete groups are descriptive: '
             'runtime censoring can depend on histories. Original model, tolerances and normalizations unchanged.')
    output.mkdir(parents=True, exist_ok=False)
    (output/'release-check-audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(result['status'], result['completed_jobs'], '/', result['expected_jobs'], flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['results', 'summary', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--archive', type=Path)
    a = parser.parse_args()
    audit(a.results.resolve(), a.summary.resolve(), a.output.resolve(), a.archive)


if __name__ == '__main__':
    main()
