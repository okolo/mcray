"""Pilot independence, geneaological aggregation and unchanged photon bins."""
import json
from pathlib import Path
import tempfile
import unittest
from run_rare_transport_pilots import jobs, begin_batch
from rare_pilot_statistics import summarize, batch_statistics
from prepare_rare_guides import products, scores


class RarePilotTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads((Path(__file__).parent/'protocols/transport_photons_v1.json').read_text())

    def test_complete_independent_pilot_matrix(self):
        plan = jobs(self.protocol, 'pilot')
        self.assertEqual(len(plan), 576)
        self.assertEqual(sum(j['histories'] for j in plan), 288000)
        seeds = {}
        for j in plan:
            self.assertGreaterEqual(j['seed'], 3200000000)
            self.assertLess(j['seed'], 2**32)
            seeds.setdefault(j['seed'], []).append(j)
        for group in seeds.values():
            self.assertEqual(len(group), 3)
            self.assertEqual(len({(j['case']['background'], j['code'], j['profile'], j['batch']) for j in group}), 1)
        self.assertTrue(set(seeds).isdisjoint({j['seed'] for j in jobs(self.protocol, 'smoke')}))

    def test_topology_is_not_a_reference_yield(self):
        inventory = [('P', 56, 26, '100000', [(55, 26), (1, 0)]),
                     ('D', 55, 26, '10000', [(55, 27)])]
        actual = scores(inventory, 'minus')
        self.assertEqual(actual[0][-1], 1.)
        self.assertEqual(actual[1][-1], 1.)
        self.assertEqual(products(11, 4, '10010', 'D'), [(1, 1), (10, 4)])
        self.assertEqual(sum(a for a, z in products(56, 26, '000001', 'P')), 56)

    def test_primary_aggregation_uses_unweighted_photon_energy_for_bins(self):
        job = dict(histories=2, seed=3210000000, profile='minus-moderate',
                   fraction=.5, copies=4, code='mcray', level=0)
        text = ('# seed=3210000000 L_mpc=1 fraction=0.5 copies=4 background_step=0.0025 eps=1e-6\n'
                '# PD_photons=On decay_photons=On primary=A56Z26;\n')
        for h in range(2):
            text += f'0 {h} A56Z26 0.5 1e20\n0 {h} A56Z26 0.5 1e20\n'
            text += f'0 {h} photon 0.125 2e14\n'
        text += '# extra_continuations=2\n# batch=0 primaries=2\n'
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'weighted.tsv'; path.write_text(text)
            result = summarize(path, job, self.protocol)
        self.assertEqual(result['means']['baryons'], 56.)
        self.assertEqual(result['means']['charge'], 26.)
        self.assertEqual(result['means']['photon_bin_2_count'], .125)
        self.assertEqual(result['means']['photon_bin_2_energy_eV'], 2e14)
        self.assertEqual(result['means']['photon_bin_1_count'], 0.)
        self.assertEqual(result['extra_continuations'], 2)
        stats = batch_statistics([dict(summary=result, seconds=1), dict(summary=result, seconds=2)])
        self.assertTrue(stats['observables']['yield:electron']['zero_variance_unresolved'])

    def test_partial_histories_are_never_recovered(self):
        job = jobs(self.protocol, 'smoke')[0]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'partial.tsv'; path.write_text('0 0 A56Z26 1 2e20\n')
            with self.assertRaises(ValueError):
                summarize(path, job, self.protocol)

    def test_started_job_without_output_is_not_relaunched(self):
        # The native adapter can fail after consuming histories but before
        # opening its final output. Absence of a TSV must not trigger a retry.
        job = jobs(self.protocol, 'smoke')[0]
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            self.assertFalse(begin_batch(out, job, ['probe', '--seed', str(job['seed'])]))
            with self.assertRaisesRegex(ValueError, 'consumed histories unknown'):
                begin_batch(out, job, ['probe', '--seed', str(job['seed'])])
            (out/(job['id']+'.tsv')).write_text('partial output')
            self.assertTrue(begin_batch(out, job, ['probe', '--seed', str(job['seed'])]))
            # Recovery eligibility does not bypass complete-history validation.
            with self.assertRaises(ValueError):
                summarize(out/(job['id']+'.tsv'), job, self.protocol)


if __name__ == '__main__':
    unittest.main()
