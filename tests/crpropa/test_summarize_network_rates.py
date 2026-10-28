#!/usr/bin/env python3
"""Ensure weak discrepancies and full/exclusive distinctions survive reporting."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class SummaryTest(unittest.TestCase):
    def test_independent_sum_weak_error_and_refinement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            rows=[dict(a=2,z=1,gamma=1e6,channel='full',independent=10.,official_direct=9.,official_gamma_interpolated=8.,refinement_absolute=1e-8),
                dict(a=2,z=1,gamma=1e6,channel='100000',independent=.05,official_direct=.04,official_gamma_interpolated=.03,refinement_absolute=1e-9),
                dict(a=5,z=2,gamma=1e6,channel='full',mcray=0.,terminal=True)]
            (root/'ref.json').write_text(json.dumps(dict(rows=rows,failed_cells=[])))
            (root/'default.tsv').write_text('2 1 1e6 full 10\n2 1 1e6 sum 0.04\n2 1 1e6 100000 0.04\n5 2 1e6 full 0\n5 2 1e6 sum 0\n')
            (root/'refined.tsv').write_text((root/'default.tsv').read_text().replace('0.04','0.05'))
            subprocess.run([sys.executable,str(Path(__file__).with_name('summarize_network_rates.py')),
                '--reference',str(root/'ref.json'),'--default',str(root/'default.tsv'),'--refined',str(root/'refined.tsv'),'--output',str(root/'out.json')],check=True,stdout=subprocess.DEVNULL)
            d=json.loads((root/'out.json').read_text());s=d['settings']['default']
            self.assertEqual(s['full']['worst_resolved_relative']['relative'],0.)
            self.assertAlmostEqual(s['sum']['worst_weak_absolute']['absolute'],-.01)
            self.assertEqual(s['sum']['worst_weak_absolute']['independent'],.05)
            self.assertEqual(d['settings']['refined']['sum']['worst_weak_absolute']['absolute'],0.)
            self.assertIn('default-to-refined',d['convergence'])


if __name__=='__main__':unittest.main()
