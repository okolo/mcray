#!/usr/bin/env python3
"""Short workflow checks without network, builds or physical calculations."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('server_runner', Path(__file__).with_name('run_tests.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RunnerTest(unittest.TestCase):
    def simulate(self, failure=None):
        with tempfile.TemporaryDirectory() as temp:
            data = Path(temp)
            fixture = data / 'input.json'
            fixture.write_text('{"seed": 42}\n')
            args = argparse.Namespace(data_dir=data, repo='https://example.invalid/repo.git',
                                      ref='talys', jobs=2, input=[fixture], workload=['approved-job'])
            calls = []

            def execute(command, cwd, logfile, env, timeout=None):
                name = logfile.stem
                calls.append(name)
                logfile.write_text(name + '\n')
                if name == failure:
                    raise subprocess.CalledProcessError(7, command)
                if failure == 'interrupted' and name == 'ctest':
                    raise KeyboardInterrupt()
                if failure == 'timeout' and name == 'ctest':
                    raise subprocess.TimeoutExpired(command, 1)
                if name == 'clone':
                    (cwd / 'source').mkdir()
                if name == 'ctest':
                    testing = cwd / 'build/Release/Testing'
                    testing.mkdir(parents=True)
                    (testing / 'LastTest.log').write_text('diagnostic')
                if name == 'workload':
                    (Path(env['MCRAY_RESULTS_DIR']) / 'spectrum.tsv').write_text('1\t2\n')

            with patch.object(runner, 'execute', side_effect=execute), patch.object(
                    runner.subprocess, 'check_output', return_value='abc123\n'), patch.object(
                    runner.platform, 'platform', return_value='test-host'):
                code = runner.run(args)
            archives = list(data.glob('*.tar.gz'))
            self.assertEqual(len(archives), 1)
            summary = json.loads(next(data.glob('2*.json')).read_text())
            self.assertEqual(summary['archive_sha256'], runner.sha256(archives[0]))
            self.assertEqual(summary['input_sha256'][str(fixture.resolve())], runner.sha256(fixture))
            with tarfile.open(archives[0]) as archive:
                members = {Path(name).name for name in archive.getnames()}
                manifest_member = next(m for m in archive.getmembers() if m.name.endswith('/manifest.json'))
                manifest = json.load(archive.extractfile(manifest_member))
            self.assertEqual(manifest['exit_code'], code)
            return code, summary, calls, members

    def test_success_freezes_commit_and_archives_workload(self):
        code, summary, calls, members = self.simulate()
        self.assertEqual(code, 0)
        self.assertEqual(summary['status'], 'passed')
        self.assertEqual(summary['commit'], 'abc123')
        self.assertIn('spectrum.tsv', members)
        self.assertIn('LastTest.log', members)
        self.assertLess(calls.index('ctest'), calls.index('workload'))

    def test_failed_build_archives_logs_and_skips_tests_and_workload(self):
        code, summary, calls, members = self.simulate('build')
        self.assertEqual(code, 7)
        self.assertEqual(summary['status'], 'failed')
        self.assertNotIn('ctest', calls)
        self.assertNotIn('workload', calls)
        self.assertIn('build.log', members)

    def test_clone_failure_still_produces_failure_archive(self):
        code, summary, calls, members = self.simulate('clone')
        self.assertEqual(code, 7)
        self.assertNotIn('commit', summary)
        self.assertEqual(calls, ['clone'])
        self.assertIn('clone.log', members)

    def test_interruption_is_not_success(self):
        code, summary, calls, members = self.simulate('interrupted')
        self.assertEqual(code, 130)
        self.assertEqual(summary['status'], 'interrupted')
        self.assertNotIn('workload', calls)

    def test_workstation_is_rejected_before_any_run(self):
        with patch.object(runner.platform, 'system', return_value='Darwin'), patch.object(
                runner.sys, 'argv', ['run_tests.py', '--server', '--data-dir', '/unused']), patch.object(
                runner, 'run') as execute_run:
            with self.assertRaises(SystemExit) as result:
                runner.main()
            self.assertEqual(result.exception.code, 2)
            execute_run.assert_not_called()

    def test_timeout_is_archived_without_starting_workload(self):
        code, summary, calls, members = self.simulate('timeout')
        self.assertEqual(code, 124)
        self.assertEqual(summary['status'], 'time_budget_exhausted')
        self.assertNotIn('workload', calls)
        self.assertIn('ctest.log', members)

    def test_timeout_terminates_descendants_and_keeps_log(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            # The child ignores TERM while its parent exits promptly. Its
            # heartbeat must stop even when the process leader has exited.
            child = directory/'child.py'
            child.write_text('import signal,time\nfrom pathlib import Path\n'
                'signal.signal(signal.SIGTERM, lambda *_: None)\n'
                'Path("child-ready").write_text("yes")\n'
                'while True:\n Path("heartbeat").write_text(str(time.monotonic()))\n time.sleep(.02)\n')
            parent = directory/'parent.py'
            parent.write_text('import subprocess,sys,time\nfrom pathlib import Path\n'
                'p=subprocess.Popen([sys.executable,"child.py"])\n'
                'while not Path("child-ready").exists(): time.sleep(.01)\n'
                'Path("child-pid").write_text(str(p.pid))\ntime.sleep(60)\n')
            with self.assertRaises(subprocess.TimeoutExpired):
                runner.execute([sys.executable, str(parent)], directory, directory/'timeout.log',
                               dict(runner.os.environ), timeout=1)
            self.assertTrue((directory/'child-ready').exists())
            heartbeat = (directory/'heartbeat').read_text()
            time.sleep(.15)
            self.assertEqual((directory/'heartbeat').read_text(), heartbeat)


if __name__ == '__main__':
    unittest.main()
