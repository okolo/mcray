"""Exercise installer failure boundaries without apt, network or a C++ build."""
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
REVISION = '34b10e27375c0e9894d7c984ac15fb96e12975ef'
TOOLS = ('prepare_talys_data.py', 'prepare_nuclear_decay_data.py',
         'prepare_gzk_data.py', 'prepare_nuclear_photon_data.py')

# These substitutes record orchestration; no scientific generator is replaced
# in the application. Real hash verification is retained in the table fixture.
PREPARER = r'''
import argparse, hashlib, json, os
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,required=True)
p.add_argument('--verify',action='store_true')
args,rest=p.parse_known_args()
out=args.output
if args.verify:
    for name in ('manifest.json','decay_manifest.json','photon_manifest.json'):
        data=json.loads((out/name).read_text())
        for file,digest in data['payload_sha256'].items():
            assert hashlib.sha256((out/file).read_bytes()).hexdigest()==digest,file
else:
    out.mkdir()
    if os.environ.get('INSTALLER_FAIL')==out.name:
        (out/'partial').write_text('incomplete')
        raise SystemExit(13)
    names=('fallback_decays.tsv','gzk_inventory.json','decays.tsv',
           'decay_stable.tsv','pd_photon_exceptions.tsv','pd_photons.tsv',
           'pd_photon_transitions.tsv','decay_photon_channels.tsv','decay_photons.tsv')
    for name in names:(out/name).write_text('fixture\n')
    payload={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in names}
    digest=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for name in ('manifest.json','decay_manifest.json','photon_manifest.json'):
        (out/name).write_text(json.dumps(dict(payload_sha256=payload,preparation_sha256=digest)))
    if out.name=='nuclear' and os.environ.get('INSTALLER_CORRUPT_STAGE'):
        (out/'decays.tsv').write_text('corrupt')
'''


class InstallerTest(unittest.TestCase):
    def setUp(self):
        # Also keep test-generated temporaries inside the project checkout.
        parent = ROOT / 'build' / 'installer-tests'
        parent.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix='path with spaces ', dir=parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.events = self.root / 'events.jsonl'
        for name in ('bin/tables', 'tools', 'fake-bin'):
            (self.root/name).mkdir(parents=True)
        shutil.copyfile(ROOT/'bin/install.sh', self.root/'bin/install.sh')
        for name in TOOLS:
            (self.root/'tools'/name).write_text(PREPARER)
        (self.root/'bin/build.sh').write_text(
            '#!/usr/bin/env bash\nexec "$(dirname "$0")/../fake-bin/build-fixture" "$@"\n')
        stub = '#!' + sys.executable + '\n' + f'''
import fcntl, json, os, shutil, sys
from pathlib import Path
name=Path(sys.argv[0]).name
args=sys.argv[1:]
root=Path({str(self.root)!r})
with (root/'events.jsonl').open('a') as f:
    f.write(json.dumps(dict(command=name,args=args,jobs=os.environ.get('MCRAY_BUILD_JOBS'),
                           tmp=os.environ.get('TMPDIR')))+'\\n')
if name=='uname':print('Linux' if args==['-s'] else 'x86_64')
elif name=='flock':
    try:fcntl.flock(9,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit(1)
elif name=='git':
    if args[0]=='clone':(Path(args[-1])/'.git').mkdir(parents=True)
    elif 'rev-parse' in args:
        print(('0'*40 if os.environ.get('INSTALLER_WRONG_REVISION') else {REVISION!r})
              if args[0]=='-C' else '1'*40)
elif name=='build-fixture' or name=='ctest':
    if os.environ.get('INSTALLER_FAIL')==name:raise SystemExit(14)
elif name=='python3' or name=='python':
    if args[:2]==['-m','venv']:
        target=Path(args[2])/'bin/python';target.parent.mkdir(parents=True)
        shutil.copyfile(sys.argv[0],target);target.chmod(0o755)
    elif args[:2]==['-m','pip']:
        if os.environ.get('INSTALLER_FAIL')=='pip':raise SystemExit(15)
        if args[2]=='freeze':print('installer-test-fixture==1')
    else:os.execv({sys.executable!r},[{sys.executable!r}]+args)
elif name=='mv':
    assert args[:2]==['-T','--no-clobber'],args
    source,dest=map(Path,args[2:])
    if os.environ.get('INSTALLER_DESTINATION_RACE'):
        dest.mkdir();(dest/'owned-by-user').write_text('preserve')
    if not dest.exists():source.rename(dest)
'''
        for name in ('uname', 'flock', 'git', 'build-fixture', 'ctest', 'python3',
                     'cmake', 'make', 'gcc', 'g++', 'gfortran', 'mv'):
            path = self.root/'fake-bin'/name
            path.write_text(stub)
            path.chmod(0o755)
        self.env = dict(os.environ, PATH=str(self.root/'fake-bin')+os.pathsep+os.environ['PATH'],
                        PYTHONDONTWRITEBYTECODE='1')
        for name in list(self.env):
            if name.startswith('INSTALLER_'):
                del self.env[name]

    def run_install(self, *args, **env):
        return subprocess.run(['bash', str(self.root/'bin/install.sh'), *args],
                              cwd=self.root, env=dict(self.env, **env),
                              text=True, capture_output=True, timeout=30)

    def events_for(self, command):
        return [x for x in map(json.loads, self.events.read_text().splitlines())
                if x['command']==command]

    @property
    def tables(self):
        return self.root/'bin/tables/nuclear'

    def test_success_preserves_other_tables_and_records_pins(self):
        sentinel = self.root/'bin/tables/sophia2'
        sentinel.write_text('existing')
        result = self.run_install('--jobs', '3')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertTrue((self.tables/'photon_manifest.json').exists())
        self.assertEqual(sentinel.read_text(), 'existing')
        self.assertEqual(self.events_for('build-fixture')[0]['jobs'], '3')
        self.assertIn(['--test-dir', str(self.root/'build/Release'),
                       '--output-on-failure', '--parallel', '3'],
                      [e['args'] for e in self.events_for('ctest')])
        git = self.events_for('git')
        self.assertTrue(any(e['args'][-2:]==['--detach', REVISION] for e in git))
        installs = [e['args'] for e in self.events_for('python') if 'install' in e['args']]
        self.assertEqual(len(installs), 1)
        self.assertIn('numpy==2.2.6', installs[0])
        self.assertIn('GitPython==3.1.45', installs[0])
        self.assertTrue(all(e['tmp']==str(self.root/'.crbeam-install/tmp')
                            for e in git if e['args'][0]!='rev-parse'))
        self.assertTrue(list((self.root/'.crbeam-install').glob('run.*/install.log')))

    def test_repeat_verifies_without_fetching_or_replacing_tables(self):
        first = self.run_install()
        self.assertEqual(first.returncode, 0, first.stdout+first.stderr)
        marker = self.tables/'decays.tsv'
        before = marker.stat().st_mtime_ns
        self.events.write_text('')
        second = self.run_install()
        self.assertEqual(second.returncode, 0, second.stdout+second.stderr)
        self.assertFalse(any(e['args'][0] in ('clone','-C') for e in self.events_for('git')))
        self.assertEqual(marker.stat().st_mtime_ns, before)
        self.assertFalse(self.events_for('python'))

    def test_invalid_existing_package_stops_before_build_and_is_preserved(self):
        self.tables.mkdir()
        (self.tables/'user-data').write_text('preserve')
        result = self.run_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.tables/'user-data').read_text(), 'preserve')
        self.assertFalse(self.events_for('build-fixture'))

    def test_corrupted_existing_package_is_not_replaced(self):
        self.assertEqual(self.run_install().returncode, 0)
        (self.tables/'decays.tsv').write_text('corrupt')
        self.events.write_text('')
        self.assertNotEqual(self.run_install().returncode, 0)
        self.assertEqual((self.tables/'decays.tsv').read_text(), 'corrupt')
        self.assertFalse(self.events_for('build-fixture'))

    def test_package_from_other_preparer_is_preserved_and_rejected(self):
        self.assertEqual(self.run_install().returncode, 0)
        manifest = self.tables/'photon_manifest.json'
        data = json.loads(manifest.read_text())
        data['preparation_sha256'] = '0'*64
        manifest.write_text(json.dumps(data))
        self.events.write_text('')
        self.assertNotEqual(self.run_install().returncode, 0)
        self.assertEqual(json.loads(manifest.read_text())['preparation_sha256'], '0'*64)
        self.assertFalse(self.events_for('build-fixture'))

    def test_existing_symlink_is_verified_without_changing_its_target(self):
        self.assertEqual(self.run_install().returncode, 0)
        immutable = self.root/'immutable-tables'
        self.tables.rename(immutable)
        self.tables.symlink_to(immutable, target_is_directory=True)
        self.events.write_text('')
        self.assertEqual(self.run_install().returncode, 0)
        self.assertTrue(self.tables.is_symlink())
        self.assertFalse(any(e['args'][0]=='clone' for e in self.events_for('git')))

    def test_second_installer_cannot_enter_the_locked_checkout(self):
        state = self.root/'.crbeam-install'
        state.mkdir()
        with (state/'install.lock').open('w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.run_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Another installer', result.stderr)
        self.assertFalse(list(state.glob('run.*')))

    def test_build_test_pip_and_preparation_failures_never_install_partial_data(self):
        for failure in ('build-fixture', 'ctest', 'pip', 'pd', 'decay', 'gzk', 'nuclear'):
            with self.subTest(failure=failure):
                result = self.run_install(INSTALLER_FAIL=failure)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.tables.exists())
        self.assertEqual(len(list((self.root/'.crbeam-install').glob('run.*'))), 7)

    def test_wrong_reference_revision_stops_before_preparation(self):
        result = self.run_install(INSTALLER_WRONG_REVISION='1')
        self.assertNotEqual(result.returncode, 0)
        logs = list((self.root/'.crbeam-install').glob('run.*/install.log'))
        self.assertIn('Wrong nuclear source revision', logs[0].read_text())
        self.assertFalse(self.tables.exists())

    def test_hash_failure_in_staged_data_never_publishes_tables(self):
        self.assertNotEqual(self.run_install(INSTALLER_CORRUPT_STAGE='1').returncode, 0)
        self.assertFalse(self.tables.exists())
        self.assertTrue(list((self.root/'.crbeam-install').glob('run.*/nuclear')))

    def test_concurrent_destination_is_preserved(self):
        result = self.run_install(INSTALLER_DESTINATION_RACE='1')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.tables/'owned-by-user').read_text(), 'preserve')
        self.assertTrue(list((self.root/'.crbeam-install').glob('run.*/nuclear')))

    def test_help_and_invalid_arguments_do_not_create_installation_state(self):
        self.assertEqual(self.run_install('--help').returncode, 0)
        for args in (('--jobs', '0'), ('--jobs',), ('--jobs', 'abc'), ('--unknown',)):
            self.assertNotEqual(self.run_install(*args).returncode, 0)
        self.assertFalse((self.root/'.crbeam-install').exists())


if __name__ == '__main__':
    unittest.main()
