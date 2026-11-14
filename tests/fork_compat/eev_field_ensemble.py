#!/usr/bin/env python3
"""Frozen EeV field/source ensemble; default fixtures never consume final seeds."""
import argparse
import concurrent.futures
import json
import math
import os
from pathlib import Path
import platform
import shlex
import signal
import subprocess
import tarfile
import threading
import time

import coupled
import source_observer_ensemble as streaming
import source_observer_pilot as pilot

HERE, ROOT = coupled.HERE, coupled.ROOT
sha256 = streaming.sha256


def design_for(mode):
    if mode not in ('fixture', 'final'):
        raise ValueError('Unknown ensemble mode')
    cases = []
    for profile, angle, theta, phi in [('pencil_z', 0, 0, 0), ('cone_z', 20, 0, 0),
                                     ('cone_rotated', 20, .4, .7)]:
        for particle, species in [('electron', 0), ('photon', 3), ('proton', 2)]:
            seed = (2010101 if mode == 'final' else 2011101)+2*len(cases)
            cases.append(dict(name=profile+'_'+particle, species=species,
                half_angle_degrees=angle, theta=theta, phi=phi, fork_seed=seed,
                modern_seed=seed+1, count=9000000 if mode == 'final' else 1000, accuracy=40))
    names = ['hit', 'arrival_0.092', 'arrival_0.094', 'arrival_0.096', 'arrival_0.098']
    names += [prefix+'_'+axis for prefix in ('position', 'final_direction', 'emitted_direction')
              for axis in 'xyz']
    return dict(schema=1, mode=mode, protocol='EEV_FIELD_ENSEMBLE.md', cases=cases,
        energy_eV=1e18, accuracy=40, alpha=.05, family_size=144,
        total_histories=sum(c['count']*2 for c in cases), observable_names=names,
        observable_widths=[1]*5+[2]*9, coordinate_scale_Mpc=.25,
        arrival_thresholds_Mpc=[.092,.094,.096,.098], source_z=.01,
        hubble_km_s_Mpc=70, omega_lambda=.7, observer_radius_Mpc=.01,
        observer_distance_Mpc=.1, duration_Mpc=.2, interactions=[], rate_path='constant')


def validate_summary(row, case, implementation):
    pilot.check_finite(row)
    n = case['count']
    if (row['case'], row['implementation'], row['histories'], row['failures'], row['emitted_weight']) != (
            case['name'], implementation, n, 0, n):
        raise ValueError('Incomplete, failed or unexpected ensemble summary')
    sums = row['sums']
    if len(sums) != 14 or sums[0] != row['hits'] or not 0 <= row['hits'] <= n:
        raise ValueError('Invalid observable counts')
    if not 0 <= sums[1] <= sums[2] <= sums[3] <= sums[4] <= sums[0]:
        raise ValueError('Invalid arrival counts')
    if any(abs(x) > n for x in sums[5:]):
        raise ValueError('Observable outside frozen bounds')
    if not 0 <= row['past_end_histories'] <= n or row['max_engine_moves'] < 1:
        raise ValueError('Invalid numerical diagnostic')


def confidence_report(rows, design):
    expected = {(c['name'], impl) for c in design['cases'] for impl in ('fork', 'modern')}
    indexed = {(r['case'], r['implementation']): r for r in rows}
    if len(rows) != len(expected) or set(indexed) != expected:
        raise ValueError('Incomplete or duplicate streams')
    for case in design['cases']:
        for impl in ('fork', 'modern'):
            validate_summary(indexed[case['name'], impl], case, impl)
    if design['mode'] == 'fixture':
        return dict(mode='fixture', streams=rows, scope='Software fixture only; no acceptance intervals.')
    if design != design_for('final'):
        raise ValueError('Frozen final design changed')
    differences, fractions = [], []
    for case in design['cases']:
        n = case['count']
        old, new = (indexed[case['name'], impl] for impl in ('fork', 'modern'))
        radius = math.sqrt(math.log(2*144/.05)/n)
        for k, (name, width) in enumerate(zip(design['observable_names'], design['observable_widths'])):
            delta, half = (new['sums'][k]-old['sums'][k])/n, width*radius
            differences.append(dict(case=case['name'], observable=name, difference=delta,
                half_width=half, interval=[delta-half, delta+half], contains_zero=abs(delta) <= half))
        for impl in ('fork', 'modern'):
            mean, half = indexed[case['name'], impl]['hits']/n, radius/math.sqrt(2)
            fractions.append(dict(case=case['name'], implementation=impl, fraction=mean,
                half_width=half, interval=[mean-half, mean+half]))
    if len(differences) != 126 or len(fractions) != 18 or differences[0]['half_width'] > .001:
        raise ValueError('Frozen family or precision changed')
    return dict(mode='final', streams=rows, confidence_level=.95, simultaneous_family_size=144,
        method='Two-sided Hoeffding plus union bound; untruncated intervals.',
        coordinate_scale_Mpc=.25, differences=differences, hit_fractions=fractions,
        detected_differences=[r for r in differences if not r['contains_zero']],
        compatible_at_reported_resolution=all(r['contains_zero'] for r in differences),
        scope='Frozen numerical algorithms/observable family only; no full physical-equivalence claim.')


def reserve(path, metadata):
    # Atomic persistent ledger: even another output directory cannot silently
    # repeat a final reservation. A failed attempt remains spent.
    with path.open('x') as ledger:
        json.dump(metadata, ledger, indent=2)
        ledger.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fork-root', 'build-dir', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--mode', choices=('fixture', 'final'), default='fixture')
    parser.add_argument('--workers', type=int, choices=range(1, 17), default=16)
    parser.add_argument('--max-seconds', type=float, default=28800)
    parser.add_argument('--ledger', type=Path)
    parser.add_argument('--command-prefix', default='', help='Prepared compiler/runtime argv prefix for fixtures')
    args = parser.parse_args()
    if not 0 < args.max_seconds <= 28800 or (args.mode == 'final' and args.ledger is None):
        parser.error('Require 0 < max-seconds <= 28800 and a persistent ledger for final collection')
    if args.mode == 'final' and args.command_prefix:
        parser.error('Final mode requires native execution; command prefixes are fixture-only')
    started = time.monotonic()
    fork, build, output = args.fork_root.resolve(), args.build_dir.resolve(), args.output.resolve()
    if output == fork or fork in output.parents:
        parser.error('Output must be outside the read-only reference')
    if output.exists() and any(output.iterdir()):
        parser.error('Preserve previous attempts; output must be new/empty')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    if args.mode == 'final' and dirty:
        parser.error('Final collection requires a clean frozen source')
    output.mkdir(parents=True, exist_ok=True); (output/'tmp').mkdir()
    design = design_for(args.mode)
    (output/'inputs.json').write_text(json.dumps(design, indent=2)+'\n')
    reference = coupled.extract(fork, output)
    fixture = output/'transport-fixture'; fixture.mkdir()
    _, grid, _ = coupled.observer.transport.inputs(fixture)
    cases = output/'ensemble.cases'
    keys = ('name', 'species', 'half_angle_degrees', 'theta', 'phi', 'fork_seed', 'modern_seed', 'count', 'accuracy')
    cases.write_text(''.join(' '.join(str(c[k]) for k in keys)+'\n' for c in design['cases']))
    protected = {p: sha256(p) for p in [output/'inputs.json', cases, grid,
        build/'mcray/libmcray.a', *[fork/n for n in reference],
        *[p for p in HERE.iterdir() if p.is_file()],
        *[p for p in (ROOT/'src/lib').iterdir() if p.suffix in ('.h', '.cpp')]]}
    metadata = dict(revision=revision, dirty=dirty, platform=platform.platform(), workers=args.workers,
        reference_sha256=reference, protected_sha256={str(p): h for p,h in protected.items()},
        max_seconds=args.max_seconds, final_reservation_spent=False)
    (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    prefix = shlex.split(args.command_prefix)
    metadata['command_prefix'] = prefix
    sources = [output/'coupled_reference.cpp', output/'reference.cpp', HERE/'eev_field_ensemble.cpp',
        HERE/'coupled_modern.cpp', *[ROOT/'src/lib'/n for n in
        ('PropagationEngine.cpp', 'Deflection3D.cpp', 'GridMagneticField.cpp', 'SpatialGeometry.cpp', 'Jet.cpp')]]
    collection = None; executor = None; timer = None; report = None
    previous_term = signal.getsignal(signal.SIGTERM)
    try:
        coupled.kernels.compile_probe(build, output, prefix, sources, 'eev-field-ensemble',
                                     timeout=min(240, args.max_seconds))
        collection = streaming.Collection(output, output/'eev-field-ensemble', cases,
                                           validator=validate_summary, prefix=prefix)
        def interrupted(signum, frame):
            collection.stop()
            raise KeyboardInterrupt('Collection interrupted')
        signal.signal(signal.SIGTERM, interrupted)
        remaining = args.max_seconds-(time.monotonic()-started)
        if remaining <= 0: raise TimeoutError('Budget exhausted before collection')
        metadata.update(executable_sha256=sha256(collection.executable),
                        start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        if args.mode == 'final':
            reserve(args.ledger.resolve(), dict(revision=revision, output=str(output), design=design,
                                              reserved_utc=metadata['start_utc']))
            metadata['final_reservation_spent'] = True
        (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
        timer = threading.Timer(remaining, collection.stop); timer.start()
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.workers)
        futures = [executor.submit(collection.run, case, impl)
                   for case in design['cases'] for impl in ('fork', 'modern')]
        rows = [f.result() for f in concurrent.futures.as_completed(futures, timeout=remaining)]
        if collection.stopped.is_set(): raise TimeoutError('Collection deadline reached')
        for path, expected in protected.items():
            if sha256(path) != expected: raise ValueError('Protected input changed: '+str(path))
        report = confidence_report(rows, design)
    except BaseException as error:
        if collection: collection.stop()
        (output/'failure.json').write_text(json.dumps(dict(error=str(error), successful_report=False), indent=2)+'\n')
        raise
    finally:
        if executor: executor.shutdown(wait=True)
        if timer: timer.cancel()
        signal.signal(signal.SIGTERM, previous_term)
        with tarfile.open(output/'metadata.tar.gz', 'w:gz') as archive:
            for path in sorted(output.iterdir()):
                if path.is_file() and path.name != 'metadata.tar.gz' and not path.name.endswith('.jsonl.gz'):
                    archive.add(path, arcname=path.name)
    elapsed = time.monotonic()-started
    if elapsed > args.max_seconds:
        raise TimeoutError('Budget exhausted before final report')
    report['elapsed_seconds'] = elapsed
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(mode=report['mode'], histories=design['total_histories'], elapsed_seconds=elapsed,
                         compatible=report.get('compatible_at_reported_resolution')), indent=2))


if __name__ == '__main__':
    main()
