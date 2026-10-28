#!/usr/bin/env python3
"""Guard full-matrix coverage and explicit failure handling in the scan driver."""
import json
from pathlib import Path
import tempfile
import unittest
from run_network_rate_scan import read_rates


class CoverageTest(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)
        (self.root/'manifest.json').write_text(json.dumps({'reachable_states':[[1,0],[1,1],[2,1],[5,2]]}))
        (self.root/'channels.tsv').write_text('2 1 100000 0 1\n')
        self.lines=[f'{a} {z} {g} {code} 0\n' for a,z,codes in [(2,1,['full','sum','100000']),(5,2,['full','sum'])] for g in [1e6,1e7] for code in codes]
        self.rates=self.root/'rates.tsv'

    def read(self,lines):
        self.rates.write_text(''.join(lines))
        return read_rates(self.rates,[1e6,1e7],self.root)

    def test_complete_including_terminals(self):
        rows,summary=self.read(self.lines)
        self.assertEqual(len(rows),10)
        self.assertEqual(summary['isotopes'],2)
        self.assertEqual(summary['missing_rows'],0)

    def test_missing_zero_row_is_not_accepted(self):
        with self.assertRaisesRegex(ValueError,'Incomplete'):self.read(self.lines[:-1])

    def test_failure_is_not_a_successful_zero(self):
        with self.assertRaisesRegex(ValueError,'Incomplete'):self.read(self.lines+['# FAILED A=2 Z=1 gamma=1e6\n'])

    def test_duplicate_negative_and_nonfinite(self):
        for lines in [self.lines+self.lines[:1], [self.lines[0].replace(' 0\n',' -1\n')]+self.lines[1:], [self.lines[0].replace(' 0\n',' nan\n')]+self.lines[1:]]:
            with self.subTest(lines=lines[:1]):
                with self.assertRaisesRegex(ValueError,'Invalid'):self.read(lines)


if __name__=='__main__':unittest.main()
