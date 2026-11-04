#!/usr/bin/env python3
"""Separate field/source cost pilots, not an acceptance test."""
import argparse
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time

import coupled
import source_observer_pilot as baseline

HERE, ROOT = coupled.HERE, coupled.ROOT
kernels = coupled.kernels


def inputs(output, mode, pilot_design="pev"):
    """Freeze the approved pilot, or a disjoint small software fixture."""
    if mode not in ('fixture', 'pilot'):
        raise ValueError('Unknown field/source mode')
    if pilot_design not in ('pev', 'eev'):
        raise ValueError('Unknown field/source design')
    count = 1000 if mode == 'pilot' else 2
    # The new EeV reservation never replays the spent PeV pilot streams.
    if pilot_design == 'eev':
        first_seed = 2007101 if mode == 'pilot' else 2008101
    else:
        first_seed = 2005101 if mode == 'pilot' else 2006101
    particles = [('electron', 0), ('photon', 3)]
    if pilot_design == 'eev':
        particles.append(('proton', 2))
    cases = []
    profiles = [('pencil_z', 0, 0, 0), ('cone_z', 20, 0, 0),
                ('cone_rotated', 20, .4, .7)]
    for profile, angle, theta, phi in profiles:
        for particle, species in particles:
            index = len(cases)
            cases.append(dict(name=profile+'_'+particle, species=species,
                              half_angle_degrees=angle, theta=theta, phi=phi,
                              fork_seed=first_seed+2*index, modern_seed=first_seed+2*index+1,
                              count=count, accuracy=20))
    path = output/'field-source.cases'
    path.write_text(''.join(' '.join(map(str, case.values()))+'\n' for case in cases))
    fixture = output/'transport-fixture'
    fixture.mkdir()
    _, grid, grid_design = coupled.observer.transport.inputs(fixture)
    design = dict(schema=1, mode=mode,
                  purpose='Cost/failure diagnostics only; no statistical equivalence verdict.',
                  cases=cases, total_histories=2*len(cases)*count, particles=[name for name,_ in particles],
                  energy_eV=1e18 if pilot_design == 'eev' else 1e15, initial_weight=1, accuracy=20,
                  cosmology=dict(hubble_km_s_Mpc=70, omega_lambda=.7, source_z=.01),
                  source_position_Mpc=[0, 0, 0], observer='strict inside sphere',
                  observer_radius_Mpc=.01, observer_distance_Mpc=.1, duration_Mpc=.2,
                  fields='Existing transport grid, zero extension',
                  grid_sha256=grid_design['sha256'][grid.name],
                  interactions=[], rate_path='constant', propagation_seed_modern=17,
                  input_sha256=kernels.digest(path), whole_run_limit_seconds=900,
                  failure_policy='Invalidate complete pilot; no retries or replacement histories.')
    (output/'inputs.json').write_text(json.dumps(design, indent=2)+'\n')
    return path, grid, design


def summarize(rows, design):
    # Reuse the complete-primary accounting and descriptive cost reducer. The
    # larger coupled state adds final directions/counters, but no confidence test.
    report = baseline.summarize(rows, design)
    for row in rows:
        if row['record'] == 'history':
            flight = row['flight']
            if row['weight'] != 1 or len(flight['direction']) != 3 or flight['moves'] < 1:
                raise ValueError('Invalid field/source weight, vector or move count')
    for group in report['groups']:
        histories = [r for r in rows if r['record'] == 'history' and
                     (r['case'], r['implementation']) == (group['case'], group['implementation'])]
        group.update(max_engine_moves=max(r['flight']['moves'] for r in histories),
                     zero_scale_calls=sum(r['flight']['zero_scale_calls'] for r in histories),
                     field_calls=sum(r['flight']['field_calls'] for r in histories),
                     seconds_per_history=group['compute_seconds']/group['histories'])
    report['mode'] = design['mode']
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root', type=Path, required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('fixture', 'pilot'), default='fixture')
    parser.add_argument('--command-prefix', default='')
    parser.add_argument('--design', choices=('pev', 'eev'), default='pev',
                        help='Original PeV pair or separately reserved EeV electron/photon/proton pilot')
    args = parser.parse_args()
    fork, build, output = args.fork_root.resolve(), args.build_dir.resolve(), args.output.resolve()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    if args.mode == 'pilot' and dirty:
        parser.error('Pilot collection requires a committed clean source revision')
    if output == fork or fork in output.parents:
        parser.error('Output must be outside the read-only reference')
    if output.exists() and any(output.iterdir()):
        parser.error('Preserve prior attempts; output must be new/empty')
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    reference = coupled.extract(fork, output)
    cases, grid, design = inputs(output, args.mode, args.design)
    prefix = shlex.split(args.command_prefix)
    metadata = dict(reference_sha256=reference,
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        dirty=dirty, host=platform.platform(), python=sys.version, command_prefix=prefix,
        adapter_sha256={p.name: kernels.digest(p) for p in HERE.iterdir() if p.is_file()},
        production_sha256={str(p.relative_to(ROOT)): kernels.digest(p)
                           for p in sorted((ROOT/'src/lib').iterdir()) if p.suffix in ('.h', '.cpp')},
        linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'))
    (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    sources = [output/'coupled_reference.cpp', output/'reference.cpp',
               HERE/'field_source_pilot.cpp', HERE/'coupled_modern.cpp',
               *[ROOT/'src/lib'/name for name in ('PropagationEngine.cpp', 'Deflection3D.cpp',
                                                'GridMagneticField.cpp', 'SpatialGeometry.cpp', 'Jet.cpp')]]
    # Compilation and execution share the approved 15-minute wall budget.
    # The executable has no child processes. Any failure leaves raw data/logs,
    # but cannot create a successful report or silently reuse reservations.
    try:
        kernels.compile_probe(build, output, prefix, sources, 'field-source-pilot', timeout=240)
        timeout = min(600, design['whole_run_limit_seconds']-(time.monotonic()-started))
        if timeout <= 0:
            raise TimeoutError('Field/source whole-run deadline expired')
        before = time.monotonic()
        with (output/'execution.log').open('w') as log, (output/'pilot.jsonl').open('w') as raw:
            subprocess.run([*prefix, str(output/'field-source-pilot'), str(cases), str(grid), args.design],
                           stdout=raw, stderr=log, check=True, timeout=timeout)
        wall = time.monotonic()-before
        rows = [json.loads(line) for line in (output/'pilot.jsonl').read_text().splitlines()]
        report = summarize(rows, design)
        report.update(process_wall_seconds=wall, whole_run_seconds=time.monotonic()-started,
                      raw_sha256=kernels.digest(output/'pilot.jsonl'),
                      executable_sha256=kernels.digest(output/'field-source-pilot'))
        if report['whole_run_seconds'] > design['whole_run_limit_seconds']:
            raise TimeoutError('Field/source whole-run deadline expired')
        (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    except Exception as error:
        (output/'failure.json').write_text(json.dumps(dict(error=str(error),
            elapsed_seconds=time.monotonic()-started, mode=args.mode), indent=2)+'\n')
        raise
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
