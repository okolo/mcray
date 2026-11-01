#!/usr/bin/env python3
"""Frozen, separate 1000-primary source/observer cost pilot; no acceptance test."""
import argparse
import json
import math
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time

import observer

HERE, ROOT = observer.HERE, observer.ROOT
kernels = observer.kernels


def inputs(output):
    # Reserve these six seeds for this pilot only. Final ensembles must use
    # different reservations; pilot counts never adapt to observed outcomes.
    profiles = [('pencil_z', 0, 0, 0), ('cone_z', 20, 0, 0),
                ('cone_rotated', 20, .4, .7)]
    cases = []
    for i, (name, angle, theta, phi) in enumerate(profiles):
        cases.append(dict(name=name, half_angle_degrees=angle, theta=theta, phi=phi,
                          fork_seed=2003101+2*i, modern_seed=2003102+2*i, count=1000))
    path = output/'pilot.cases'
    path.write_text(''.join(' '.join(map(str, case.values()))+'\n' for case in cases))
    design = dict(schema=1, purpose='Cost and failure diagnostics only; no acceptance verdict.',
                  cases=cases, total_histories=6000,
                  cosmology=dict(hubble_km_s_Mpc=70, omega_lambda=.7, source_z=.01),
                  particle='photon', energy_eV=1e12, initial_weight=1,
                  source_position_Mpc=[0, 0, 0], observer='strict inside sphere',
                  observer_radius_Mpc=.01, observer_distance_Mpc=.1,
                  observer_axis='Declared source axis; common represented center in both paths.',
                  duration_Mpc=.2, rate_path='constant', interactions=[], fields=False,
                  propagation_seed_modern=17, # No interaction draws in this model.
                  failure_policy='Invalidate the entire pilot; no replacement histories.',
                  input_sha256=kernels.digest(path))
    (output/'inputs.json').write_text(json.dumps(design, indent=2)+'\n')
    return path, design


def check_finite(value):
    if isinstance(value, dict):
        for item in value.values():
            check_finite(item)
    elif isinstance(value, list):
        for item in value:
            check_finite(item)
    elif isinstance(value, (float, int)) and not math.isfinite(value):
        raise ValueError('Nonfinite pilot output')


def summarize(rows, design):
    """Require every primary exactly once before deriving descriptive totals."""
    expected = {(case['name'], implementation, i)
                for case in design['cases'] for implementation in ('fork', 'modern')
                for i in range(case['count'])}
    seen, groups, timings = set(), {}, {}
    for row in rows:
        check_finite(row)
        group = row['case'], row['implementation']
        if row['record'] == 'timing':
            if group in timings or row['compute_seconds'] < 0:
                raise ValueError('Duplicate or invalid pilot timing')
            timings[group] = row['compute_seconds']
            continue
        if row['record'] != 'history':
            raise ValueError('Unknown pilot record')
        key = (*group, row['primary'])
        if key in seen or key not in expected:
            raise ValueError('Duplicate or unexpected pilot primary')
        if type(row['flight']['hit']) is not bool or type(row['flight']['past_end']) is not bool:
            raise ValueError('Invalid pilot hit/time flags')
        if len(row['direction']) != 3 or len(row['flight']['position_Mpc']) != 3:
            raise ValueError('Invalid pilot vector')
        seen.add(key)
        groups.setdefault(group, []).append(row)
    if seen != expected or set(timings) != set(groups):
        raise ValueError('Incomplete pilot; partial or failed histories are not misses')
    summaries = []
    for (case, implementation), histories in groups.items():
        hits = [row for row in histories if row['flight']['hit']]
        summaries.append(dict(case=case, implementation=implementation, histories=len(histories),
            hits=len(hits), misses=len(histories)-len(hits),
            emitted_weight=math.fsum(row['weight'] for row in histories),
            hit_weight=math.fsum(row['weight'] for row in hits),
            zero_weight_histories=sum(row['weight'] == 0 for row in histories),
            negative_weight_histories=sum(row['weight'] < 0 for row in histories),
            max_direction_norm_error=max(abs(math.sqrt(math.fsum(x*x for x in row['direction']))-1)
                                         for row in histories),
            past_end_histories=sum(row['flight']['past_end'] for row in histories),
            max_positive_end_error_Mpc=max(0, max(row['flight']['end_error_Mpc'] for row in histories)),
            max_historical_steps=max(row['historical_steps'] for row in histories),
            compute_seconds=timings[case, implementation]))
    return dict(histories=len(seen), groups=summaries,
                failures=0, # A failed engine/process/partial result cannot reach this point.
                scope=design['purpose'],
                timing_scope='Sampling plus propagation; excludes JSON/compile/startup. '
                             'Historical propagation uses a reduced interaction-free driver.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root', type=Path, required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True,
                        help='New, empty directory; protects spent pilot reservations and logs')
    parser.add_argument('--command-prefix', default='')
    args = parser.parse_args()
    fork, output, build = args.fork_root.resolve(), args.output.resolve(), args.build_dir.resolve()
    if output == fork or fork in output.parents:
        parser.error('Output must be outside the reference')
    if output.exists() and any(output.iterdir()):
        parser.error('Output must be empty; preserve earlier pilot attempts')
    output.mkdir(parents=True, exist_ok=True)
    reference = observer.extract(fork, output)
    cases, design = inputs(output)
    prefix = shlex.split(args.command_prefix)
    metadata = dict(reference_sha256=reference,
        modern_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        modern_working_tree_dirty=bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True)),
        adapter_sha256={p.name: kernels.digest(p) for p in sorted(HERE.iterdir()) if p.is_file()},
        production_sha256={str(p.relative_to(ROOT)): kernels.digest(p)
                           for p in sorted((ROOT/'src/lib').iterdir()) if p.suffix in ('.h', '.cpp')},
        linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'),
        environment=dict(python=sys.version, host=platform.platform(), command_prefix=prefix))
    (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    sources = [output/'observer_reference.cpp', output/'reference.cpp',
               HERE/'source_observer_pilot.cpp', HERE/'observer_modern.cpp',
               *[ROOT/'src/lib'/name for name in ('GridMagneticField.cpp', 'Deflection3D.cpp',
                                                'SpatialGeometry.cpp', 'PropagationEngine.cpp')]]
    kernels.compile_probe(build, output, prefix, sources, 'source-observer-pilot')
    before = time.monotonic()
    # Keep partial output for diagnosis, but publish a report only on complete success.
    with (output/'execution.log').open('w') as log, (output/'pilot.jsonl').open('w') as raw:
        subprocess.run([*prefix, str(output/'source-observer-pilot'), str(cases)],
                       stdout=raw, stderr=log, check=True, timeout=180)
    wall = time.monotonic()-before
    rows = [json.loads(line) for line in (output/'pilot.jsonl').read_text().splitlines()]
    report = summarize(rows, design)
    report.update(process_wall_seconds=wall,
                  executable_sha256=kernels.digest(output/'source-observer-pilot'),
                  raw_sha256=kernels.digest(output/'pilot.jsonl'))
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
