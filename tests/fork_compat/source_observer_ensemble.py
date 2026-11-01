#!/usr/bin/env python3
"""Frozen source/observer comparison with bounded memory and compressed raw histories."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import tarfile
import threading
import time

import observer

HERE, ROOT = observer.HERE, observer.ROOT
kernels = observer.kernels


def sha256(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def design_for(mode):
    design = json.loads((HERE/'source_observer_ensemble.json').read_text())
    if mode == 'fixture':
        # Replay only the already spent pilot streams. Final reservations are
        # never consumed by a smoke test; no acceptance report for fixtures.
        for i, case in enumerate(design['cases']):
            case.update(count=1000, fork_seed=2003101+2*i, modern_seed=2003102+2*i)
        design['total_histories'] = 6000
    design['mode'] = mode
    return design


def validate_summary(row, case, implementation):
    def finite(value):
        if isinstance(value, dict):
            for item in value.values(): finite(item)
        elif isinstance(value, list):
            for item in value: finite(item)
        elif isinstance(value, (int, float)) and not math.isfinite(value):
            raise ValueError('Nonfinite ensemble summary')
    finite(row)
    n = case['count']
    if (row['case'], row['implementation'], row['histories'], row['failures'], row['emitted_weight']) != (
            case['name'], implementation, n, 0, n):
        raise ValueError('Incomplete, failed or unexpected ensemble summary')
    sums = row['sums']
    if len(sums) != 11 or sums[0] != row['hits'] or not 0 <= row['hits'] <= n:
        raise ValueError('Invalid observable counts')
    if not 0 <= sums[1] <= sums[2] <= sums[3] <= sums[4] <= sums[0]:
        raise ValueError('Invalid arrival counts')
    if any(abs(x) > n for x in sums[5:]):
        raise ValueError('Observable outside frozen bounds')
    if not 0 <= row['past_end_histories'] <= n:
        raise ValueError('Invalid diagnostic counts')


def confidence_report(rows, design):
    """Compute the fixed simultaneous family only after all six streams pass."""
    expected = {(case['name'], impl) for case in design['cases'] for impl in ('fork', 'modern')}
    by_key = {(row['case'], row['implementation']): row for row in rows}
    if len(rows) != len(expected) or set(by_key) != expected:
        raise ValueError('Incomplete or duplicate streams')
    for case in design['cases']:
        for impl in ('fork', 'modern'):
            validate_summary(by_key[case['name'], impl], case, impl)
    if design['mode'] != 'final':
        return dict(mode='fixture', streams=rows, scope='Fixture replay only; no acceptance intervals.')
    differences, fractions = [], []
    for case in design['cases']:
        old, modern = (by_key[case['name'], impl] for impl in ('fork', 'modern'))
        n = case['count']
        radius = math.sqrt(math.log(2*design['family_size']/design['alpha'])/n)
        for k, (name, width) in enumerate(zip(design['observable_names'], design['observable_widths'])):
            delta = modern['sums'][k]/n-old['sums'][k]/n
            half = width*radius
            differences.append(dict(case=case['name'], observable=name, difference=delta,
                half_width=half, interval=[delta-half, delta+half], contains_zero=abs(delta) <= half))
        for impl in ('fork', 'modern'):
            mean = by_key[case['name'], impl]['hits']/n
            half = radius/math.sqrt(2)
            fractions.append(dict(case=case['name'], implementation=impl, fraction=mean,
                                  half_width=half, interval=[mean-half, mean+half]))
    if len(differences)+len(fractions) != 39:
        raise ValueError('Frozen confidence family changed')
    if any(row['half_width'] > .001 for row in differences if row['observable'] == 'hit'):
        raise ValueError('Insufficient hit-fraction precision')
    return dict(mode='final', streams=rows, confidence_level=.95, simultaneous_family_size=39,
        method='Two-sided Hoeffding bounds with union bound; untruncated intervals.',
        differences=differences, hit_fractions=fractions,
        detected_differences=[row for row in differences if not row['contains_zero']],
        compatible_at_reported_resolution=all(row['contains_zero'] for row in differences),
        scope='Frozen observables only; no exact or full physical-equivalence claim.')


class Collection:
    """Each child owns one RNG stream; no shared historical iostream/RNG state."""
    def __init__(self, output, executable, cases):
        self.output, self.executable, self.cases = output, executable, cases
        self.lock, self.processes, self.stopped = threading.Lock(), [], threading.Event()

    def stop(self):
        with self.lock:
            self.stopped.set()
            for process in self.processes:
                if process.poll() is None:
                    try: os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError: pass

    def run(self, case, implementation):
        stem = self.output/(case['name']+'-'+implementation)
        raw_path = stem.with_suffix('.jsonl.gz')
        summary_path = stem.with_suffix('.json')
        raw_hash, lines = hashlib.sha256(), 0
        env = dict(os.environ, TMPDIR=str(self.output/'tmp'))
        with stem.with_suffix('.log').open('wb') as log:
            with self.lock:
                if self.stopped.is_set(): raise RuntimeError('Collection cancelled')
                process = subprocess.Popen([str(self.executable), str(self.cases), case['name'],
                    implementation, str(summary_path)], stdout=subprocess.PIPE, stderr=log,
                    cwd=ROOT, env=env, start_new_session=True)
                self.processes.append(process)
            try:
                with raw_path.open('wb') as disk, gzip.GzipFile(fileobj=disk, mode='wb',
                        filename='', mtime=0, compresslevel=1) as archive:
                    while True:
                        block = process.stdout.read(1024*1024)
                        if not block: break
                        if self.stopped.is_set(): raise TimeoutError('Collection deadline reached')
                        raw_hash.update(block); lines += block.count(b'\n')
                        archive.write(block)
                if process.wait() != 0: raise RuntimeError('Failed stream: '+stem.name)
            finally:
                process.stdout.close()
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL); process.wait()
        if lines != case['count']: raise ValueError('Missing or excess raw primary records')
        row = json.loads(summary_path.read_text())
        validate_summary(row, case, implementation)
        # Read back the compressed artifact, verifying its CRC, complete bytes
        # and hash without loading a history collection into memory.
        verified = hashlib.sha256()
        with gzip.open(raw_path, 'rb') as archive:
            for block in iter(lambda: archive.read(1024*1024), b''):
                if self.stopped.is_set(): raise TimeoutError('Archive verification deadline reached')
                verified.update(block)
        if verified.digest() != raw_hash.digest(): raise ValueError('Compressed history verification failed')
        row.update(raw_sha256=raw_hash.hexdigest(), gzip_sha256=sha256(raw_path),
                   gzip_bytes=raw_path.stat().st_size, raw_lines=lines)
        stem.with_suffix('.verified.json').write_text(json.dumps(row, indent=2)+'\n')
        print(stem.name+': complete and archive verified', flush=True)
        return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root', type=Path, required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('fixture', 'final'), default='fixture')
    parser.add_argument('--workers', type=int, choices=range(1, 7), default=6)
    parser.add_argument('--max-seconds', type=int, default=14400)
    args = parser.parse_args()
    started = time.monotonic()
    if not 0 < args.max_seconds <= 14400: parser.error('Limit must be 1..14400 seconds')
    fork, build, output = args.fork_root.resolve(), args.build_dir.resolve(), args.output.resolve()
    if output == fork or fork in output.parents: parser.error('Output must be outside the reference')
    if output.exists() and any(output.iterdir()): parser.error('Preserve previous attempt; output must be empty')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True))
    if args.mode == 'final' and dirty: parser.error('Final collection requires a clean frozen revision')
    output.mkdir(parents=True, exist_ok=True); (output/'tmp').mkdir()
    design = design_for(args.mode)
    (output/'inputs.json').write_text(json.dumps(design, indent=2)+'\n')
    reference = observer.extract(fork, output)
    cases = output/'ensemble.cases'
    cases.write_text(''.join(' '.join(map(str, [c[key] for key in
        ('name', 'half_angle_degrees', 'theta', 'phi', 'fork_seed', 'modern_seed', 'count')]))+'\n'
        for c in design['cases']))
    metadata = dict(revision=revision, dirty=dirty, platform=platform.platform(), workers=args.workers,
        reference_sha256=reference, input_sha256=sha256(output/'inputs.json'),
        adapter_sha256={p.name: sha256(p) for p in sorted(HERE.iterdir()) if p.is_file()},
        production_sha256={str(p.relative_to(ROOT)): sha256(p) for p in sorted((ROOT/'src/lib').iterdir())
                           if p.suffix in ('.cpp', '.h')},
        linked_core_sha256=sha256(build/'mcray/libmcray.a'),
        local_limit_seconds=args.max_seconds, final_reservation_spent=False)
    (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    sources = [output/'observer_reference.cpp', output/'reference.cpp', HERE/'source_observer_ensemble.cpp',
        HERE/'observer_modern.cpp', *[ROOT/'src/lib'/name for name in
            ('GridMagneticField.cpp', 'Deflection3D.cpp', 'SpatialGeometry.cpp', 'PropagationEngine.cpp')]]
    kernels.compile_probe(build, output, [], sources, 'source-observer-ensemble',
                          timeout=max(1, args.max_seconds-(time.monotonic()-started)))
    collection = Collection(output, output/'source-observer-ensemble', cases)
    remaining = args.max_seconds-(time.monotonic()-started)
    if remaining <= 0: raise TimeoutError('Time budget exhausted before collection')
    timer = threading.Timer(remaining, collection.stop)
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.workers)
    try:
        metadata.update(final_reservation_spent=args.mode == 'final',
                        executable_sha256=sha256(collection.executable),
                        start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        (output/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
        timer.start()
        futures = [executor.submit(collection.run, case, impl)
                   for case in design['cases'] for impl in ('fork', 'modern')]
        rows = [future.result() for future in concurrent.futures.as_completed(futures, timeout=remaining)]
        if collection.stopped.is_set(): raise TimeoutError('Collection deadline reached')
        for name, expected in metadata['adapter_sha256'].items():
            if sha256(HERE/name) != expected: raise ValueError('Adapter changed during collection')
        for name, expected in metadata['production_sha256'].items():
            if sha256(ROOT/name) != expected: raise ValueError('Production source changed during collection')
        if sha256(build/'mcray/libmcray.a') != metadata['linked_core_sha256']:
            raise ValueError('Linked core changed during collection')
        for name, expected in reference.items():
            if sha256(fork/name) != expected: raise ValueError('Reference changed during collection')
        report = confidence_report(rows, design)
        report.update(elapsed_seconds=time.monotonic()-started)
        (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    except BaseException as error:
        collection.stop()
        (output/'failure.json').write_text(json.dumps(dict(error=str(error),
            elapsed_seconds=time.monotonic()-started, successful_report=False), indent=2)+'\n')
        raise
    finally:
        timer.cancel(); executor.shutdown(wait=True)
        with tarfile.open(output/'metadata.tar.gz', 'w:gz') as archive:
            for path in sorted(output.iterdir()):
                if path.is_file() and path.name != 'metadata.tar.gz' and not path.name.endswith('.jsonl.gz'):
                    archive.add(path, arcname=path.name)
    print(json.dumps(dict(mode=report['mode'], elapsed_seconds=report['elapsed_seconds'],
        compatible=report.get('compatible_at_reported_resolution'),
        detected_differences=len(report.get('detected_differences', []))), indent=2))


if __name__ == '__main__':
    main()
