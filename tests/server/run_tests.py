#!/usr/bin/env python3
"""Fetch a frozen remote revision and archive a server test run, including failures."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tarfile
import time
import uuid


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def execute(command, cwd, logfile, env, timeout=None):
    with logfile.open('wb') as stream:
        # A workload can spawn worker processes. Terminate its whole process
        # group on budget expiry so abandoned transport jobs cannot continue.
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=5)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            # The leader can exit before a descendant that ignored TERM.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            raise
        if code:
            raise subprocess.CalledProcessError(code, command)


def archive_run(run_dir, data_dir, manifest):
    results = run_dir / 'results'
    manifest['finished_utc'] = utc_now()
    manifest['payload_sha256'] = {
        str(p.relative_to(results)): sha256(p)
        for p in sorted(results.rglob('*')) if p.is_file()
    }
    (results / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    archive = data_dir / (run_dir.name + '.tar.gz')
    partial = archive.with_suffix(archive.suffix + '.partial')
    with tarfile.open(partial, 'w:gz') as output:
        output.add(results, arcname=run_dir.name)
    partial.replace(archive)
    summary = dict(manifest, archive=archive.name, archive_sha256=sha256(archive))
    (data_dir / (run_dir.name + '.json')).write_text(json.dumps(summary, indent=2) + '\n')
    return archive


def run(args):
    started = time.monotonic()
    maximum = getattr(args, 'max_runtime_seconds', None)
    data_dir = args.data_dir.expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    run_dir = data_dir / 'work' / run_id
    results = run_dir / 'results'
    results.mkdir(parents=True)
    checkout = run_dir / 'source'
    env = dict(os.environ, MCRAY_RESULTS_DIR=str(results),
               MCRAY_BUILD_JOBS=str(args.jobs))
    manifest = {'started_utc': utc_now(), 'status': 'running', 'repo': args.repo,
                'requested_ref': args.ref, 'host': platform.platform(),
                'architecture': platform.machine(), 'jobs': args.jobs, 'steps': [],
                'runner_sha256': sha256(Path(__file__)),
                'workload': args.workload,
                'input_sha256': {str(p.resolve()): sha256(p) for p in args.input}}
    if maximum is not None:
        manifest['max_runtime_seconds'] = maximum
    code = 0
    try:
        commands = [
            ('clone', ['git', 'clone', '--no-checkout', '--', args.repo, str(checkout)], run_dir),
            ('fetch', ['git', 'fetch', 'origin', args.ref], checkout),
            ('checkout', ['git', 'checkout', '--detach', 'FETCH_HEAD'], checkout),
            ('build', ['bash', 'bin/build.sh', 'Release'], checkout),
            ('ctest', ['ctest', '--test-dir', 'build/Release', '--output-on-failure',
                       '--output-junit', str(results / 'ctest.xml')], checkout),
        ]
        if args.workload:
            commands.append(('workload', args.workload, checkout))
        for name, command, cwd in commands:
            step = {'name': name, 'command': command, 'started_utc': utc_now()}
            manifest['steps'].append(step)
            print(f'{run_id}: {name}', flush=True)
            remaining = None if maximum is None else maximum-(time.monotonic()-started)
            if remaining is not None and remaining <= 0:
                raise subprocess.TimeoutExpired(command, maximum)
            execute(command, cwd, results / (name + '.log'), env, timeout=remaining)
            step['exit_code'] = 0
            step['finished_utc'] = utc_now()
            if name == 'checkout':
                manifest['commit'] = subprocess.check_output(
                    ['git', 'rev-parse', 'HEAD'], cwd=checkout, text=True).strip()
                execute(['git', 'show', '--no-patch', '--format=fuller', 'HEAD'],
                        checkout, results / 'source-revision.txt', env)
        manifest['status'] = 'passed'
    except subprocess.TimeoutExpired as error:
        code = 124
        manifest.update(status='time_budget_exhausted', error=str(error))
        manifest['steps'][-1]['exit_code'] = code
    except subprocess.CalledProcessError as error:
        code = error.returncode if error.returncode > 0 else 128 - error.returncode
        manifest['status'] = 'failed'
        manifest['error'] = str(error)
        manifest['steps'][-1]['exit_code'] = code
    except KeyboardInterrupt:
        code = 130
        manifest['status'] = 'interrupted'
    except Exception as error:
        code = 1
        manifest['status'] = 'failed'
        manifest['error'] = str(error)
    finally:
        testing = checkout / 'build/Release/Testing'
        if testing.exists():
            shutil.copytree(testing, results / 'ctest-diagnostics', dirs_exist_ok=True)
        manifest['exit_code'] = code
        archive = archive_run(run_dir, data_dir, manifest)
        print(f"{manifest['status']}: {archive}", flush=True)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--server', action='store_true',
                        help='Confirm execution on an external server, not the workstation')
    parser.add_argument('--repo', default='https://github.com/okolo/mcray.git')
    parser.add_argument('--ref', default='talys', help='Remote branch, tag or commit (default talys)')
    parser.add_argument('--data-dir', required=True, type=Path)
    parser.add_argument('--jobs', default=2, type=int)
    parser.add_argument('--max-runtime-seconds', type=int,
                        help='Cap clone/build/test/workload time; always archive afterwards')
    parser.add_argument('--input', action='append', type=Path, default=[],
                        help='Input file to hash; repeat for every external input/manifest')
    parser.add_argument('workload', nargs=argparse.REMAINDER,
                        help='Optional command after --; cwd is the frozen checkout')
    args = parser.parse_args()
    if not args.server or platform.system() != 'Linux' or platform.machine() != 'x86_64':
        parser.error('Requires --server and an external native Linux amd64 server')
    if args.jobs < 1 or args.ref.startswith('-'):
        parser.error('Require positive --jobs and a ref that does not start with -')
    if args.max_runtime_seconds is not None and args.max_runtime_seconds < 1:
        parser.error('--max-runtime-seconds must be positive')
    if args.workload[:1] == ['--']:
        args.workload = args.workload[1:]
    for path in args.input:
        if not path.is_file():
            parser.error(f'Input is not a file: {path}')
    return run(args)


if __name__ == '__main__':
    sys.exit(main())
