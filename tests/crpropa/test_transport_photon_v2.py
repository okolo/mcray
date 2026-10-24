"""Preparation guards: preserve the physical family and account for all work."""
import json
from pathlib import Path
import unittest
from prepare_transport_photon_v2 import ordinary_plan, RARE


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads((Path(__file__).parent/'protocols/transport_photons_v1.json').read_text())
        self.times = {(bg, code, level): {'seconds_per_primary': .01*(level+1)*(100 if code == 'mcray' else 1)}
                      for bg in self.protocol['backgrounds']
                      for code in ['mcray', 'crpropa'] for level in range(3)}

    def test_budget_and_family(self):
        before = json.dumps(self.protocol, sort_keys=True)
        plan = ordinary_plan(self.protocol, self.times)
        self.assertEqual(json.dumps(self.protocol, sort_keys=True), before)
        self.assertEqual(set(plan['ordinary_observables']) | set(RARE),
                         set(self.protocol['targets']['cmb/mcray']['means']))
        self.assertFalse(set(plan['ordinary_observables']) & set(RARE))
        self.assertLessEqual(plan['expected_ordinary_wall_hours']*1.15, 20.)
        sizes, seeds = {}, {}
        for job in plan['jobs']:
            group = job['case']['background'], job['level'], job['code']
            sizes[group] = sizes.get(group, 0)+job['histories']
            seeds.setdefault(job['seed'], []).append(job)
            self.assertGreaterEqual(job['seed'], 3000000000)
            self.assertLess(job['seed'], 2**32)
        for bg in self.protocol['backgrounds']:
            for level in range(3):
                self.assertEqual(sizes[bg, level, 'mcray'], sizes[bg, level, 'crpropa'])
        # Seed reuse is allowed only for numerical pairing of the same primary batch.
        for entries in seeds.values():
            self.assertEqual(len({(j['case']['background'], j['code'], j['batch']) for j in entries}), 1)
            self.assertEqual(len({j['level'] for j in entries}), len(entries))

    def test_incomplete_or_unusable_timing_matrix(self):
        missing = dict(self.times)
        missing.pop(next(iter(missing)))
        with self.assertRaises(ValueError):
            ordinary_plan(self.protocol, missing)
        with self.assertRaises(ValueError):
            ordinary_plan(self.protocol, self.times, hours=1e-10)
        with self.assertRaises(ValueError):
            ordinary_plan(self.protocol, self.times, workers=0)


if __name__ == '__main__':
    unittest.main()
