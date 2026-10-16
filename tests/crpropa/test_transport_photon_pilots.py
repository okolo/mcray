#!/usr/bin/env python3
"""Verify ensemble identity, history bookkeeping and failure reporting."""
import json
from pathlib import Path
import tempfile
import unittest
from run_transport_photon_pilots import make_jobs,summarize,reports,sha


class TransportPilotsTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'particles.tsv'
        self.job=make_jobs('initial','smoke')[0];self.job['histories']=2

    def write(self,body,footer=True):
        self.path.write_text(body+('# batch=0 primaries=2\n' if footer else ''))

    def test_matrix_and_paired_seeds(self):
        jobs=make_jobs('initial','pilot');self.assertEqual(len(jobs),1440)
        self.assertEqual(len({j['id'] for j in jobs}),1440)
        seeds={}
        for j in jobs:seeds.setdefault((j['case']['id'],j['code'],j['batch']),set()).add(j['seed'])
        self.assertTrue(all(len(v)==1 for v in seeds.values()))
        self.assertEqual(len(set().union(*seeds.values())),len(seeds))
        self.assertTrue(all(j['histories']==1000 for j in jobs))
        self.assertEqual(len(make_jobs('full','pilot')),8640)

    def test_complete_correlated_history(self):
        self.write('0 0 A56Z26 1 56e18\n0 1 A55Z26 1 55e18\n0 1 neutron 1 1e18\n')
        d=summarize(self.path,self.job)
        self.assertEqual(d['survival'],.5);self.assertEqual(d['baryons'],56);self.assertEqual(d['charge'],26)

    def test_photons_and_spectra(self):
        self.job['case']['pd_photons']=True
        self.write('0 0 A56Z26 1 56e18\n0 0 photon 1 1e15\n0 0 photon 1 2e15\n0 1 A56Z26 1 56e18\n')
        d=summarize(self.path,self.job);self.assertEqual(d['yield:photon'],1);self.assertEqual(d['photon_histories'],.5)
        self.assertEqual(d['energy:photon'],1.5e15)

    def test_reject_incomplete_wrong_ids_or_bookkeeping(self):
        bodies=['0 0 A56Z26 1 56e18\n','0 0 A56Z26 1 56e18\n0 2 A56Z26 1 56e18\n',
            '0 0 A55Z26 1 55e18\n0 1 A56Z26 1 56e18\n','0 0 A56Z26 1 nan\n0 1 A56Z26 1 56e18\n',
            '0 0 A56Z26 1 56e18\n0 0 photon 1 1e15\n0 1 A56Z26 1 56e18\n']
        for body in bodies:
            self.write(body)
            with self.subTest(body=body),self.assertRaises(ValueError):summarize(self.path,self.job)
        self.write('0 0 A56Z26 1 56e18\n0 1 A56Z26 1 56e18\n',False)
        with self.assertRaises(ValueError):summarize(self.path,self.job)

    def test_empty_report_remains_incomplete(self):
        reports(self.root,dict(jobs=[self.job]))
        d=json.loads((self.root/'pilot-report.json').read_text())
        self.assertEqual(d['completed_jobs'],0);self.assertEqual(d['expected_jobs'],1)
        self.assertIn('no agreement verdict',d['status'])

    def test_paired_refinement_and_independent_program_comparison(self):
        all_jobs=make_jobs('initial','pilot');case=all_jobs[0]['case']['id']
        jobs=[j for j in all_jobs if j['case']['id']==case and j['batch']<2 and j['level']<2]
        for job in jobs:
            path=self.root/(job['id']+'.tsv');path.write_text('retained output\n')
            mean=2*job['batch']+job['level']+(3 if job['code']=='crpropa' else 0)
            (self.root/(job['id']+'.done.json')).write_text(json.dumps(dict(output_sha256=sha(path),observables={'yield:photon':mean})))
        reports(self.root,dict(jobs=jobs));d=json.loads((self.root/'pilot-report.json').read_text())
        self.assertEqual(d['completed_jobs'],8)
        self.assertTrue(all(r['paired_difference']==1 and r['paired_standard_error']==0 for r in d['numerical_refinements']))
        self.assertTrue(all(r['difference']==-3 and r['combined_standard_error']>0 for r in d['program_comparisons']))


if __name__=='__main__':unittest.main()
