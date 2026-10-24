"""Analytic weight identities and finite-path, competing-hazard controls."""
import math
import statistics
import unittest
from rare_channel_sampling import (HistoryRewards, competing_decay_fixture,
                                   draw_channel, proposal, split_weight)


class RareChannelTests(unittest.TestCase):
    def test_two_transition_path_expectation(self):
        # Exhaustively integrate the changed branch choices, not MC noise.
        p = [.99999, .00001]
        q = proposal(p, [1., 100000.], .8)
        mean = 0.
        for i in range(2):
            for j in range(2):
                mid_i = math.fsum(q[:i])+q[i]/2
                mid_j = math.fsum(q[:j])+q[j]/2
                a = draw_channel(p, q, mid_i)
                b = draw_channel(p, q, mid_j)
                mean += q[i]*q[j]*a.likelihood*b.likelihood*(i+j)
        self.assertAlmostEqual(mean, 2*p[1], places=16)

    def test_no_self_normalization_and_genealogy(self):
        rewards = HistoryRewards()
        for w in split_weight(.6, 3):
            rewards.add('count', w, 2.)
        self.assertAlmostEqual(rewards.values()['count'], 1.2)

    def test_finite_path_with_competition(self):
        # Expected unconditional yield includes both hazards and no-event paths.
        rows = competing_decay_fixture(2., 3., .2, [.99999, .00001],
                                       [1., 100000.], .8, 4, 16000, 341281)
        values = [r.get('1', 0.) for r in rows]
        exact = 2/5*(-math.expm1(-1.))*.00001
        se = statistics.stdev(values)/math.sqrt(len(values))
        self.assertLess(abs(statistics.mean(values)-exact), 6*se)
        # Same model, different guide; expected value must remain unchanged.
        rows = competing_decay_fixture(2., 3., .2, [.99999, .00001],
                                       [1., 100000.], .5, 2, 16000, 341282)
        values = [r.get('1', 0.) for r in rows]
        se = statistics.stdev(values)/math.sqrt(len(values))
        self.assertLess(abs(statistics.mean(values)-exact), 6*se)

    def test_zero_rates_and_full_support(self):
        self.assertEqual(competing_decay_fixture(0., 0., 1., [1., 0.],
                         [1., 10.], .8, 2, 3, 12), [{}, {}, {}])
        self.assertEqual(proposal([1., 0.], [0., 10.], .8), (1., 0.))
        with self.assertRaises(ValueError):
            draw_channel([.9, .1], [1., 0.], .5)
        with self.assertRaises(ValueError):
            proposal([.9, .1], [1., 2.], 1.)


if __name__ == '__main__':
    unittest.main()
