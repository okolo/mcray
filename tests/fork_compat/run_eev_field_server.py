#!/usr/bin/env python3
"""Single-entry Ubuntu 22.04 amd64 EeV comparison launcher (eight hours total)."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import threading
import time
import uuid

REFERENCE = 'e1fc71f613efe7fecf3210e15339c34814245cd4'
GTEST_SHA256 = '6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5'
LIMIT = 8*3600


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024*1024), b''): digest.update(block)
    return digest.hexdigest()


def tree_pids(pid):
    """Include workers that deliberately started their own process sessions."""
    result = [pid]
    try:
        children = Path('/proc/{}/task/{}/children'.format(pid, pid)).read_text().split()
    except FileNotFoundError:
        children = []
    for child in children: result.extend(tree_pids(int(child)))
    return result


def kill_tree(process, captured=None):
    # Capture descendants before killing the leader, which would reparent them.
    if captured is None and process.poll() is not None: return
    for pid in reversed(captured if captured is not None else tree_pids(process.pid)):
        try: os.kill(pid, signal.SIGKILL)
        except ProcessLookupError: pass
    process.wait()


class Supervisor:
    def __init__(self, directory, seconds=LIMIT):
        self.directory, self.started = directory, time.monotonic()
        self.seconds, self.active, self.lock = seconds, [], threading.RLock()

    def remaining(self, reserve=None):
        if reserve is None: reserve = min(120, self.seconds/10)
        value = self.seconds-(time.monotonic()-self.started)-reserve
        if value <= 0: raise TimeoutError('Eight-hour budget exhausted')
        return value

    def stop(self):
        with self.lock:
            for process in self.active: kill_tree(process)

    def execute(self, name, command, cwd, env):
        timeout = self.remaining()
        print(name, flush=True)
        with (self.directory/(name+'.log')).open('wb') as log:
            with self.lock:
                process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                self.active.append(process)
            try:
                code = process.wait(timeout=timeout)
            except BaseException:
                # Give the ensemble its TERM handler briefly, then ensure all
                # descendants (including separate worker sessions) are gone.
                captured = tree_pids(process.pid)
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=5)
                except (ProcessLookupError, subprocess.TimeoutExpired): pass
                kill_tree(process, captured)
                raise
            if code: raise subprocess.CalledProcessError(code, command)


def extract_archive(archive, destination):
    """Extract only regular files/dirs; reject traversal and archive links."""
    destination.mkdir()
    with tarfile.open(archive, 'r:gz') as source:
        members = source.getmembers()
        for member in members:
            path = (destination/member.name).resolve()
            if destination.resolve() not in path.parents and path != destination.resolve():
                raise ValueError('Unsafe archive path')
            if not (member.isdir() or member.isfile()): raise ValueError('Archive links are not accepted')
        source.extractall(destination, members=members)
    folders = list(destination.iterdir())
    if len(folders) != 1 or not folders[0].is_dir(): raise ValueError('Expected one archive root')
    return folders[0]


def run(args):
    data = args.data_dir.expanduser().resolve(); data.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8]
    directory = data/run_id; directory.mkdir(); (directory/'tmp').mkdir()
    supervisor = Supervisor(directory)
    env = dict(os.environ, TMPDIR=str(directory/'tmp'), PYTHONDONTWRITEBYTECODE='1',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), requested_revision=args.revision,
        repository=args.repository, runner_sha256=sha256(Path(__file__)), host=platform.platform(),
        workers=args.workers, fixture_only=args.fixture_only, max_seconds=LIMIT, status='running')
    ledger = data/'eev-field-ensemble-v1.reservation.json'
    def save():
        temporary = directory/'server-status.pending.json'
        temporary.write_text(json.dumps(manifest, indent=2)+'\n')
        temporary.replace(directory/'server-status.json')
    def hard_deadline():
        supervisor.stop()
        manifest.update(status='hard_deadline', exit_code=124)
        save()
        os._exit(124)
    def interrupted(signum, frame):
        supervisor.stop()
        raise KeyboardInterrupt('Server launcher interrupted')
    previous_signals = {s: signal.signal(s, interrupted) for s in (signal.SIGTERM, signal.SIGINT)}
    timer = threading.Timer(LIMIT, hard_deadline); timer.daemon = True; timer.start()
    code = 1
    try:
        save()
        if not args.fixture_only and ledger.exists():
            raise FileExistsError('Final reservation already spent; preserve ledger and previous results')
        if not args.fixture_only and shutil.disk_usage(data).free < 64*2**30:
            raise RuntimeError('Final run requires at least 64 GiB free in the test-data directory')
        checkout = directory/'source'
        supervisor.execute('clone', ['git', 'clone', '--no-checkout', '--', args.repository, str(checkout)], directory, env)
        # A bundle can already contain this object without advertising it as a ref.
        available = subprocess.run(['git', 'cat-file', '-e', args.revision+'^{commit}'], cwd=checkout,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
        if not available:
            supervisor.execute('fetch', ['git', 'fetch', 'origin', args.revision], checkout, env)
        supervisor.execute('checkout', ['git', 'checkout', '--detach', args.revision], checkout, env)
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=checkout, text=True).strip()
        if commit != args.revision: raise ValueError('Resolved revision differs from request')
        candidate_launcher = checkout/'tests/fork_compat/run_eev_field_server.py'
        if sha256(candidate_launcher) != manifest['runner_sha256']:
            raise ValueError('Launcher differs from the selected commit; transfer the paired script and REVISION')
        manifest['commit'] = commit; save()
        supervisor.execute('versions', ['bash', '-c',
            'uname -a; cat /etc/os-release; g++ --version; gfortran --version; cmake --version; '
            'python3 --version; dpkg-query -W libgsl-dev libboost-filesystem-dev libboost-system-dev libyaml-cpp-dev libssl-dev'], checkout, env)
        def download(name, url, supplied):
            if supplied: return supplied.expanduser().resolve()
            path = directory/(name+'.tar.gz')
            supervisor.execute(name+'-download', ['curl', '--fail', '--location', '--retry', '3',
                '--connect-timeout', '30', '--output', str(path), url], directory, env)
            return path
        archive = download('reference', 'https://codeload.github.com/okolo/mcray/tar.gz/'+REFERENCE, args.reference_archive)
        manifest['reference_archive_sha256'] = sha256(archive)
        reference = extract_archive(archive, directory/'reference')
        gtest_archive = download('googletest', 'https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.18.0', args.gtest_archive)
        if sha256(gtest_archive) != GTEST_SHA256: raise ValueError('GoogleTest archive hash mismatch')
        gtest = extract_archive(gtest_archive, directory/'googletest')
        build = checkout/'build/server-eev'
        supervisor.execute('configure', ['cmake', '-S', 'src/app/crbeam', '-B', str(build),
            '-DCMAKE_BUILD_TYPE=Release', '-DBUILD_TESTING=ON', '-DCMAKE_C_COMPILER=gcc',
            '-DCMAKE_CXX_COMPILER=g++', '-DCMAKE_Fortran_COMPILER=gfortran',
            '-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST='+str(gtest)], checkout, env)
        supervisor.execute('build', ['cmake', '--build', str(build), '--parallel', str(args.workers)], checkout, env)
        supervisor.execute('ctest', ['ctest', '--test-dir', str(build), '--output-on-failure',
            '--parallel', str(min(4,args.workers)), '--output-junit', str(directory/'ctest.xml')], checkout, env)
        supervisor.execute('python-guards', ['python3', 'tests/test_fork_transport.py', '-v'], checkout, env)
        common = ['python3', 'tests/fork_compat/eev_field_ensemble.py', '--fork-root', str(reference),
                  '--build-dir', str(build), '--workers', str(args.workers)]
        supervisor.execute('fixture', [*common, '--mode', 'fixture', '--output', str(directory/'fixture'),
            '--max-seconds', str(min(600, supervisor.remaining()))], checkout, env)
        if not args.fixture_only:
            supervisor.execute('ensemble', [*common, '--mode', 'final', '--output', str(directory/'ensemble'),
                '--ledger', str(ledger), '--max-seconds', str(supervisor.remaining())], checkout, env)
            result = json.loads((directory/'ensemble/report.json').read_text())
            manifest['compatible_at_reported_resolution'] = result['compatible_at_reported_resolution']
        manifest['reference_revision'] = REFERENCE
        manifest['status'] = 'fixture_passed' if args.fixture_only else 'completed'
        code = 0
    except BaseException as error:
        supervisor.stop()
        code = 124 if isinstance(error, (TimeoutError, subprocess.TimeoutExpired)) else 1
        manifest.update(status='failed', error=str(error))
    finally:
        manifest.update(exit_code=code, elapsed_seconds=time.monotonic()-supervisor.started,
                        finished_utc=datetime.now(timezone.utc).isoformat())
        save()
        # Raw histories already have streaming gzip/verified hashes. Archive
        # only small server metadata, avoiding a second copy of large payloads.
        try:
            with tarfile.open(directory/'server-metadata.tar.gz', 'w:gz') as archive:
                for path in sorted(directory.iterdir()):
                    if path.is_file() and path.name != 'server-metadata.tar.gz': archive.add(path, arcname=path.name)
        except Exception as error:
            code = 1
            manifest.update(status='failed', exit_code=code, error='Metadata archive: '+str(error))
            save()
        finally:
            timer.cancel()
            for signum, handler in previous_signals.items(): signal.signal(signum, handler)
    print(json.dumps(dict(directory=str(directory), **manifest), indent=2))
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--server', action='store_true', help='Confirm external-server execution')
    parser.add_argument('--repository', default='https://github.com/okolo/mcray.git', help='Git URL or transferred Git bundle')
    parser.add_argument('--revision', required=True, help='Exact published or bundled 40-character commit')
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, choices=range(1,17), default=16)
    parser.add_argument('--fixture-only', action='store_true', help='Full setup/checks, but never consume final seeds')
    parser.add_argument('--reference-archive', type=Path, help='Optional transferred immutable snapshot tar.gz')
    parser.add_argument('--gtest-archive', type=Path, help='Optional SHA-verified GoogleTest 1.18.0 tar.gz')
    args = parser.parse_args()
    release = Path('/etc/os-release')
    ubuntu = release.is_file() and 'ID=ubuntu' in release.read_text() and 'VERSION_ID="22.04"' in release.read_text()
    if not args.server or platform.system() != 'Linux' or platform.machine() != 'x86_64' or not ubuntu:
        parser.error('Requires --server on Ubuntu 22.04 amd64 (fixture-only may run in a Linux test container)')
    if not re.fullmatch('[0-9a-f]{40}', args.revision): parser.error('Require the exact lowercase 40-character commit')
    sys.exit(run(args))


if __name__ == '__main__': main()
