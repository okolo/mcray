#!/usr/bin/env python3
"""Guard acceptance bookkeeping, history statistics and decision boundaries."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from run_transport_photon_acceptance import PROTOCOL, jobs, summary, interval, report, reserve_campaign, sha, command


class AcceptanceTest(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads(PROTOCOL.read_text())
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_cap_and_paired_refinement_seeds(self):
        all_jobs = jobs(self.protocol)
        self.assertEqual(len(all_jobs), 272)
        self.assertEqual(len({j['id'] for j in all_jobs}), 272)
        for bg in ['cmb', 'cmb-inoue']:
            for code in ['mcray', 'crpropa']:
                subset = [j for j in all_jobs if j['code'] == code and j['case']['background'] == bg]
                self.assertEqual(sum(j['histories'] for j in subset), 136000)
                for batch in range(10):
                    self.assertEqual(len({j['seed'] for j in subset if j['batch'] == batch}), 1)
        independent = {(j['code'], j['case']['background'], j['batch']): j['seed'] for j in all_jobs}
        self.assertEqual(len(set(independent.values())), len(independent))
        altered = copy.deepcopy(self.protocol); altered['main_batches'] += 1
        with self.assertRaisesRegex(ValueError, 'cap'):
            jobs(altered)
        self.assertTrue(set(independent.values()).isdisjoint(
            set(range(2026101600, 2026101600+1200)) | set(range(2027101600, 2027101600+1200))))

    def test_fixed_photon_edges_and_correlated_descendants(self):
        job = jobs(self.protocol)[0]; job['histories'] = 2
        path = self.root/'rows.tsv'
        path.write_text(f'# process=talys-decay L_mpc=1 seed={job["seed"]} accuracy=1e-6 background_step_decades=.0025\n'
                        '# PD_photons=On decay_photons=On\n'
                        '0 0 A56Z26 1 56e18\n0 0 photon 1 1e14\n0 0 photon 1 1e15\n'
                        '0 1 A55Z26 1 55e18\n0 1 neutron 1 1e18\n# batch=0 primaries=2\n')
        result = summary(path, job, self.protocol)
        self.assertEqual(result['survival'], .5)
        self.assertEqual(result['photon_bin_0_count'], 0.)
        self.assertEqual(result['photon_bin_1_count'], .5)
        self.assertEqual(result['photon_bin_2_count'], .5)
        self.assertEqual(result['energy:bound'], 55.5e18)
        self.assertEqual(result['yield:A55Z26'], .5)
        path.write_text('0 0 A56Z26 1 56e18\n')
        with self.assertRaises(ValueError):
            summary(path, job, self.protocol)

    def test_equivalence_disagreement_and_zero_samples(self):
        def check(values, target=1., epsilon=0.):
            return interval(values, target, 1., epsilon, False, self.protocol, 384)
        self.assertEqual(check([1.]*48)['status'], 'inconclusive')
        self.assertEqual(check([0.]*48)['status'], 'inconclusive')
        self.assertEqual(check([.9999, 1.0001]*24)['status'], 'agreement')
        self.assertEqual(check([1.02, 1.0201]*24)['status'], 'disagreement')
        self.assertEqual(check([.9999, 1.0001]*24, epsilon=.011)['status'], 'inconclusive')
        self.assertEqual(interval([.9999, 1.0001]*24, 1., 0., 0., False,
                                  self.protocol, 384)['status'], 'inconclusive')

    def test_duplicate_campaign_is_rejected_without_overwriting_reservation(self):
        ledger = self.root/'ledger.json'
        reservation = dict(identity={'protocol': 'fixed'}, output='/first', histories=136000)
        reserve_campaign(ledger, reservation)
        before = ledger.read_bytes()
        reserve_campaign(ledger, reservation)
        changed = dict(reservation, output='/second')
        with self.assertRaisesRegex(ValueError, 'already reserved'):
            reserve_campaign(ledger, changed)
        self.assertEqual(ledger.read_bytes(), before)

    def test_full_report_recomputes_raw_histories_and_rejects_changed_statistics(self):
        protocol = copy.deepcopy(self.protocol)
        protocol.update(main_batches=2, refinement_batches=2, histories_per_batch=2)
        all_jobs = jobs(protocol)
        for job in all_jobs:
            c = protocol['native_controls'][job['level']]
            step = protocol['reference_step_fractions'][job['level']]
            header = (f'# process=talys-decay L_mpc=1 seed={job["seed"]} accuracy={c["eps"]} '
                      f'background_step_decades={c["step"]} step_fraction={step}\n'
                      '# PD_photons=On decay_photons=On\n')
            rows = ('0 0 A56Z26 1 56e18\n0 1 A56Z26 1 56e18\n' if job['batch'] == 0 else
                    '0 0 A56Z26 1 56e18\n0 1 A55Z26 1 55e18\n0 1 neutron 1 1e18\n')
            target = self.root/(job['id']+'.tsv')
            target.write_text(header+rows+'# batch=0 primaries=2\n')
            marker = dict(job=job, output_sha256=sha(target), observables=summary(target, job, protocol))
            (self.root/(job['id']+'.done.json')).write_text(json.dumps(marker))
        result = report(self.root, protocol, all_jobs)
        self.assertTrue(result['execution_complete'])
        self.assertEqual(len(result['observables']), 384)
        self.assertNotEqual(result['status'], 'agreement')
        marker_path = self.root/(all_jobs[0]['id']+'.done.json')
        marker = json.loads(marker_path.read_text()); marker['observables']['survival'] += .1
        marker_path.write_text(json.dumps(marker))
        with self.assertRaisesRegex(ValueError, 'raw histories'):
            report(self.root, protocol, all_jobs)

    def test_missing_batches_never_accept(self):
        result = report(self.root, self.protocol, jobs(self.protocol))
        self.assertEqual(result['status'], 'inconclusive')
        self.assertFalse(result['execution_complete'])
        self.assertEqual(result['simultaneous_family_size'], 384)

    def test_reference_command_uses_the_retained_runtime_library(self):
        bundle = self.root/'bundle'
        (bundle/'nuclei').mkdir(parents=True)
        (bundle/'nuclei/masses.tsv').write_text('56 26 52100 NIST\n')
        references = self.root/'relocated-reference'
        job = next(j for j in jobs(self.protocol) if j['code'] == 'crpropa')
        argv, env = command(job, self.protocol, bundle, references, None,
                            references/'reference-build/crpropa_decay_probe', self.root/'out.tsv')
        self.assertEqual(env['LD_LIBRARY_PATH'].split(':')[0],
                         str(references/'reference-build/reference'))
        self.assertEqual(env['CRPROPA_DATA_PATH'], str(references/'cmb'))
        self.assertEqual(argv[0], str(references/'reference-build/crpropa_decay_probe'))


if __name__ == '__main__':
    unittest.main()
