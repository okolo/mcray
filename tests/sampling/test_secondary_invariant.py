import csv
import io
from pathlib import Path
import tempfile
import unittest

from invariant_diagnostics import PROBABILITIES, summarize


class InvariantControls(unittest.TestCase):
    def test_bounded_enclosure_and_tail_controls_are_required_together(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ICS.tsv';self.fixture(path)
            text=path.read_text()+'# envelope upper=1.005 lower=1 cells=8 evaluations=10 checks=100 max_ratio=0.999 max_local_ratio=0.999\n'
            text+=''.join(f'# benchmark bounded {r} 128 0.015 0\n' for r in range(3))
            path.write_text(text)
            self.assertEqual(summarize([path])['status'],'unresolved')
            path.write_text(text+'# bounded_tail_proposals 19\n')
            self.assertEqual(summarize([path])['status'],'within_control')
            path.write_text((text+'# bounded_tail_proposals 19\n').replace('upper=1.005','upper=0.5'))
            self.assertEqual(summarize([path])['status'],'unresolved')

    def fixture(self, path, error=0., reference_change=0.):
        out=io.StringIO()
        out.write('# total=2 total_error_estimate=1e-10 reference=2 reference_coarse=2\n')
        writer=csv.writer(out,delimiter='\t')
        writer.writerow(('u','x','small_tail_reference','relative_probability_error','reference_change_relative','evaluations','iterations'))
        for i,u in enumerate(PROBABILITIES):
            writer.writerow((u,i+1,min(u,1-u),error,reference_change,10,2))
        for repeat in range(3):
            for mode,seconds in [('cdf',.03),('rejection',.01),('nr',.02)]:
                out.write(f'# benchmark {mode} {repeat} 128 {seconds} 0\n')
        path.write_text(out.getvalue())

    def test_complete_controls_preserve_full_construction_cost(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ICS.tsv';self.fixture(path)
            result=summarize([path]);self.assertEqual(result['status'],'within_control')
            self.assertAlmostEqual(result['cases']['ICS']['cdf_over_rejection_time'],3)

    def test_reference_uncertainty_cannot_be_reported_as_accuracy(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ICS.tsv';self.fixture(path,reference_change=2e-4)
            result=summarize([path]);self.assertEqual(result['status'],'unresolved')
            self.assertEqual(result['cases']['ICS']['status'],'inconclusive_reference')

    def test_material_quantile_error_fails_even_with_valid_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ICS.tsv';self.fixture(path,error=.01)
            self.assertEqual(summarize([path])['cases']['ICS']['status'],'deviation')

    def test_nonfinite_and_incomplete_diagnostics_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ICS.tsv';self.fixture(path,error=float('nan'))
            with self.assertRaises(ValueError):summarize([path])
            self.fixture(path);path.write_text(path.read_text().replace('# benchmark cdf 2 128 0.03 0\n',''))
            with self.assertRaises(ValueError):summarize([path])


if __name__=='__main__':unittest.main()
