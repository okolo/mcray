"""Release limits must not turn uncertain observations into accepted physics."""
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from run_release_check import make_plan, execute, finish_report, PROTOCOL, RARE


class ReleaseCheckTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads(PROTOCOL.read_text())
        self.times = {(bg, code, level): {'seconds_per_primary': .01*(level+1)*(100 if code == 'mcray' else 1)}
                      for bg in self.protocol['backgrounds'] for code in ['mcray', 'crpropa'] for level in range(3)}

    def test_new_fixed_campaign_keeps_physics_and_reserves_time(self):
        before = json.dumps(self.protocol, sort_keys=True)
        plan = make_plan(self.protocol, self.times, 16)
        self.assertEqual(before, json.dumps(self.protocol, sort_keys=True))
        self.assertEqual(len(plan['jobs']), 1024)
        self.assertLessEqual(plan['expected_ordinary_wall_hours']*1.25+.5, 8)
        self.assertEqual(plan['failure_probability_primary'], .008)
        seeds = {}
        for j in plan['jobs']:
            self.assertGreaterEqual(j['seed'], 3300000000)
            self.assertLess(j['seed'], 2**32)
            seeds.setdefault(j['seed'], []).append(j)
        for entries in seeds.values():
            self.assertEqual(len({(j['case']['background'], j['code'], j['batch']) for j in entries}), 1)
            self.assertEqual(len(entries), len({j['level'] for j in entries}))

    def test_deadline_never_starts_or_replays_a_job(self):
        job = make_plan(self.protocol, self.times, 16)['jobs'][0]
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            with patch('run_release_check.command', return_value=(['unused'], {})), \
                    patch('run_release_check.subprocess.run') as run:
                self.assertFalse(execute(job, self.protocol, out, out, out, out, out, time.time()-1))
                run.assert_not_called()
                self.assertFalse(list(out.iterdir()))
                (out/(job['id']+'.tsv')).write_text('untrusted existing output')
                with self.assertRaises(ValueError):
                    execute(job, self.protocol, out, out, out, out, out, time.time()+60)
                run.assert_not_called()

    def test_rare_and_partial_evidence_cannot_establish_full_acceptance(self):
        rare = dict(observable=RARE[0], status='agreement')
        normal = dict(observable='survival', status='agreement')
        fake = dict(status='agreement', observables=[rare, normal], execution_complete=False,
                    note='', main_family_status='agreement')
        with tempfile.TemporaryDirectory() as directory, patch('run_release_check.report', return_value=fake):
            r = finish_report(Path(directory), self.protocol, [])
            self.assertEqual(r['ordinary_comparison_status'], 'inconclusive')
            self.assertFalse(r['full_physical_acceptance'])
            self.assertEqual(rare['status'], 'inconclusive')
            self.assertEqual(rare['diagnostic_status'], 'agreement')


if __name__ == '__main__':
    unittest.main()
