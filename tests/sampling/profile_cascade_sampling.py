#!/usr/bin/env python3
"""Bounded transport timing, with source, executable, commands and outputs frozen."""
import argparse
import hashlib
import json
import os
import pathlib
import platform
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    parser.add_argument('--primaries', type=int, default=30)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--scope', choices=('cascade', 'interactions'), default='cascade')
    parser.add_argument('--distances', type=float, nargs='+', default=[.01, .1])
    parser.add_argument('--rate-variations', type=float, nargs='+', default=[0.])
    args = parser.parse_args()
    if not (0 < args.seconds <= 1800 and args.primaries > 0 and args.repeats >= 3):
        parser.error('Require 0 < seconds <= 1800, positive primaries and repeats >= 3')
    if any(not 0 < d <= 1 for d in args.distances) or any(not 0 <= v <= .1 for v in args.rate_variations):
        parser.error('Distances must be in (0,1] Mpc and rate variations in [0,0.1]')
    if args.scope == 'cascade':
        scenarios = [(f'd{d}-v{v}', [str(d), str(v)]) for d in args.distances for v in args.rate_variations]
    else:
        scenarios = [(f'{p}-{e:g}', [p, str(e)]) for p in ('ICS', 'PP')
                     for e in ([1e13, 1e14, 1e15, 1e16, 1e18, 1e20] if p == 'ICS' else [1e14, 1e15, 1e16, 1e18, 1e20])]
    binary = args.binary.resolve()
    repo = pathlib.Path(__file__).resolve().parents[2]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    deadline = time.monotonic() + args.seconds
    records = []
    manifest = dict(status='running', platform=platform.platform(), machine=platform.machine(),
                    head=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    arguments=vars(args).copy(), records=records,
                    scope=args.scope+'; single-thread loop only; excludes initialization and final output. '
                          'Invariant selection includes its normalization. Same corrected secondary laws in both modes.')
    manifest['arguments'] = {k: str(v) if isinstance(v, pathlib.Path) else v for k, v in manifest['arguments'].items()}
    (output / 'source.patch').write_bytes(subprocess.check_output(['git', '-C', str(repo), 'diff', 'HEAD']))
    subprocess.run(['git', '-C', str(repo), 'archive', '-o', str(output / 'source.tar'), 'HEAD'], check=True)
    (output / 'probe.cpp').write_bytes((repo / 'tests/sampling/probe.cpp').read_bytes())
    (output / 'profile_cascade_sampling.py').write_bytes(pathlib.Path(__file__).read_bytes())
    (output / 'em_sampling_probe').write_bytes(binary.read_bytes())
    (output / 'em_sampling_probe').chmod(0o755)
    try:
        for scenario, tail in scenarios:
            for repeat in range(args.repeats):
                # Alternate order to reduce systematic warm-up and thermal effects.
                for mode in (('rejection', 'current') if repeat % 2 == 0 else ('current', 'rejection')):
                    stem = f'{scenario}-r{repeat}-{mode}'
                    result = output / (stem + '.tsv')
                    command = [str(binary), args.scope, mode, str(args.primaries), str(20261020 + repeat),
                               str(result), *tail]
                    with (output / (stem + '.log')).open('w') as log:
                        subprocess.run(command, env=dict(os.environ, MCRAY_PROFILE_S='1'), stdout=log,
                                       stderr=subprocess.STDOUT, check=True, timeout=max(.001, min(60., deadline-time.monotonic())))
                    timing = next(line for line in result.read_text().splitlines() if line.startswith('# transport_seconds='))
                    fields = dict(item.split('=') for item in timing[2:].split())
                    total = float(fields['transport_seconds'])
                    invariant = float(fields['invariant_seconds'])
                    records.append(dict(scenario=scenario, repeat=repeat,
                                        mode=mode, command=command, transport_seconds=total,
                                        invariant_seconds=invariant, invariant_fraction=invariant/total,
                                        invariant_calls=int(fields['invariant_calls'])))
                    (output / 'timings.json').write_text(json.dumps(manifest, indent=2)+'\n')
                    print(stem, f'total={total:.4f}s S={invariant/total:.2%}', flush=True)
        summaries = []
        for scenario, _ in scenarios:
            groups = {mode: [r for r in records if r['mode'] == mode and r['scenario'] == scenario] for mode in ('rejection', 'current')}
            medians = {mode: {key: statistics.median(r[key] for r in group)
                             for key in ('transport_seconds', 'invariant_seconds', 'invariant_fraction', 'invariant_calls')}
                       for mode, group in groups.items()}
            summaries.append(dict(scenario=scenario, medians=medians,
                                  transport_ratio=medians['current']['transport_seconds']/medians['rejection']['transport_seconds']))
        manifest.update(status='complete', summaries=summaries)
    except BaseException as error:
        manifest.update(status='execution_error', error=str(error))
        raise
    finally:
        (output / 'timings.json').write_text(json.dumps(manifest, indent=2)+'\n')
        checksums = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir()) if p.is_file()}
        (output / 'sha256.json').write_text(json.dumps(checksums, indent=2)+'\n')


if __name__ == '__main__':
    main()
