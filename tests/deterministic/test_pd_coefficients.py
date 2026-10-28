#!/usr/bin/env python3
"""Runtime representation rules that must not be silently normalized away."""
import math
import unittest
from pd_coefficients import products,nearest,effective_branches,ratio_probability


class CoefficientTest(unittest.TestCase):
    def test_full_rate_not_exclusive_sum_is_event_rate(self):
        from transport import solve
        model=dict(schema='mcray-deterministic-transport-v1',initial='He4',energy_scale_eV=4.,states=[
            dict(id='He4',a=4,z=2,energy_eV=4.,rate_per_Mpc=6.,channels=[dict(probability=2/3,products={'D':2}),dict(probability=1/3,products={'D':2})]),
            dict(id='D',a=2,z=1,energy_eV=2.,rate_per_Mpc=0.,channels=[])])
        self.assertAlmostEqual(solve(model,1.)['initial_state_yield'],math.exp(-6.),13)

    def test_products_and_multiple_alphas(self):
        residual,emitted=products((56,26),'000002')
        self.assertEqual(residual,(48,22));self.assertEqual(emitted[4,2],2)
        self.assertEqual(sum(a*v for (a,z),v in emitted.items()),56)
        with self.assertRaises(ValueError):products((2,1),'900000')

    def test_nearest_node_uses_cpp_half_up(self):
        self.assertEqual(nearest(10**6.5),13)
        self.assertEqual(nearest(1e6),0);self.assertEqual(nearest(1e14),200)
        with self.assertRaises(ValueError):nearest(1e15)

    def test_serialized_branch_tail_follows_runtime(self):
        self.assertEqual(effective_branches([.2,.3,.499999999]),[.2,.3,.5])
        self.assertEqual(effective_branches([.7,.4,.1]),[.7,.30000000000000004,0.])
        with self.assertRaises(ValueError):effective_branches([-.1,1.1])

    def test_saturation_and_zero_denominator(self):
        self.assertEqual(ratio_probability(2.,1.),1.)
        self.assertEqual(ratio_probability(1.,0.),1.)
        self.assertEqual(ratio_probability(0.,0.),0.)
        with self.assertRaises(ValueError):ratio_probability(-1.,2.)


if __name__=='__main__':unittest.main()
