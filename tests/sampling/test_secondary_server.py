import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from run_secondary_server import run, seed
class ServerBoundaries(unittest.TestCase):
    def test_deadline_terminates_process_and_retains_log(self):
        with tempfile.TemporaryDirectory() as d:
            log=Path(d)/'log';begin=time.monotonic()
            with self.assertRaises(TimeoutError):
                run([sys.executable,'-u','-c','import time; print("preserved"); time.sleep(60)'],log,begin+.2)
            self.assertLess(time.monotonic()-begin,3)
            self.assertIn('preserved',log.read_text())
    def test_failed_command_retains_stderr(self):
        with tempfile.TemporaryDirectory() as d:
            log=Path(d)/'log'
            with self.assertRaises(RuntimeError):
                run([sys.executable,'-c','import sys; print("diagnostic",file=sys.stderr); sys.exit(4)'],log,time.monotonic()+3)
            self.assertIn('diagnostic',log.read_text())
    def test_seeds_are_fixed_distinct_and_nonzero(self):
        labels=[str(i)+mode for i in range(1000) for mode in ('old','current','pilot')]
        values=[seed(label) for label in labels]
        self.assertEqual(values,[seed(label) for label in labels])
        self.assertEqual(len(values),len(set(values)))
        self.assertTrue(all(0<v<2**32 for v in values))

class FailureRetention(unittest.TestCase):
    def test_preflight_failure_keeps_status_and_attempts_archive(self):
        import contextlib
        import io
        import json
        from unittest.mock import patch
        import run_secondary_server as runner
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);nr=root/'nr';nr.mkdir();output=root/'output'
            argv=['runner','--server','--repo','unused','--nr',str(nr),'--data-dir',str(output),'--hours','.01']
            commands=[]
            def fake_run(command,log,deadline,cwd=None):
                commands.append(command)
                if command[0]=='tar':Path(command[command.index('-cf')+1]).write_bytes(b'archive')
            with patch.object(sys,'argv',argv),patch.object(runner.platform,'system',return_value='Linux'),patch.object(runner.platform,'machine',return_value='x86_64'),patch.object(runner.shutil,'which',return_value='/tool'),patch.object(runner,'run',side_effect=fake_run),contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(runner.main(),1)
            directory=next(p for p in output.iterdir() if p.is_dir())
            final=json.loads((directory/'final-status.json').read_text())
            self.assertEqual(final['status'],'execution_error')
            self.assertIn('nr3.h',final['error'])
            self.assertTrue((directory/'sha256.json').is_file())
            self.assertTrue(directory.with_suffix('.tar.gz').is_file())
            self.assertTrue(directory.with_suffix('.status.json').is_file())
            self.assertEqual(commands[0][0],'tar')


class NRProvisioning(unittest.TestCase):
    def fixture(self, directory):
        import hashlib
        import json
        import subprocess
        from run_secondary_server import NR_HEADERS
        root=Path(directory);repo=root/'upstream';repo.mkdir();tools=root/'tools';tools.mkdir()
        def git(*args):
            return subprocess.check_output(['git','-C',str(repo),'-c','user.name=Fixture','-c','user.email=fixture@localhost','-c','commit.gpgsign=false','-c','core.hooksPath=/dev/null',*args],stderr=subprocess.PIPE,text=True).strip()
        git('init','-q')
        hashes={}
        for name in NR_HEADERS:
            path=repo/name;path.write_text('historical '+name+'\n')
            hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        git('add','.')
        git('commit','-qm','test: create NR fixture')
        commit=git('rev-parse','HEAD')
        (tools/'nr_reference_hashes.json').write_text(json.dumps({'commit':commit,'headers':hashes}))
        return root,repo,tools,git,commit
    def test_download_uses_pinned_commit_not_branch_head(self):
        import json
        from run_secondary_server import prepare_nr
        with tempfile.TemporaryDirectory() as d:
            root,repo,tools,git,commit=self.fixture(d)
            (repo/'nr3.h').write_text('changed branch head\n')
            git('commit','-qam','test: advance moving NR branch')
            result=root/'run';result.mkdir();(result/'logs').mkdir()
            headers=prepare_nr(None,result,tools,time.monotonic()+10,str(repo))
            self.assertEqual((headers/'nr3.h').read_text(),'historical nr3.h\n')
            provenance=json.loads((result/'nr-provenance.json').read_text())
            self.assertEqual(provenance['resolved_commit'],commit)
            self.assertTrue((result/'logs/nr-fetch.log').is_file())
    def test_local_override_never_downloads(self):
        from unittest.mock import patch
        import run_secondary_server as runner
        with tempfile.TemporaryDirectory() as d:
            root,repo,tools,git,commit=self.fixture(d)
            result=root/'run';result.mkdir()
            with patch.object(runner,'run') as command:
                runner.prepare_nr(repo,result,tools,time.monotonic()+5)
                command.assert_not_called()
            self.assertFalse((result/'nr-upstream').exists())
    def test_hash_mismatch_fails_before_build_and_preserves_provenance(self):
        from run_secondary_server import prepare_nr
        with tempfile.TemporaryDirectory() as d:
            root,repo,tools,git,commit=self.fixture(d)
            (repo/'nr3.h').write_text('wrong header\n')
            result=root/'run';result.mkdir()
            with self.assertRaisesRegex(RuntimeError,'hash mismatch: nr3.h'):
                prepare_nr(repo,result,tools,time.monotonic()+5)
            self.assertTrue((result/'nr-provenance.json').is_file())
            self.assertEqual((result/'nr-private/nr3.h').read_text(),'wrong header\n')

if __name__=='__main__':unittest.main()
