#!/usr/bin/env python3
"""Analytic controls for independently propagated means, not acceptance budgets."""
import copy
import math
import unittest
import numpy as np
from transport import solve,line_moments


def fixture(rate=2.):
    return dict(schema='mcray-deterministic-transport-v1',initial='D',energy_scale_eV=100.,photon_bin_edges_eV=[0,2,4,None],states=[
        dict(id='D',a=2,z=1,energy_eV=100.,rate_per_Mpc=rate,channels=[dict(probability=1.,products={'p':1,'n':1},photons=[dict(mean_energy_eV=2.,probability=.25)])]),
        dict(id='p',a=1,z=1,energy_eV=50.,rate_per_Mpc=0.,channels=[]),dict(id='n',a=1,z=0,energy_eV=50.,rate_per_Mpc=0.,channels=[])])


class TransportTest(unittest.TestCase):
    def test_exponential_survival_and_integrated_rewards(self):
        d=fixture();r=solve(d,.7);survival=math.exp(-1.4)
        self.assertAlmostEqual(r['initial_state_yield'],survival,13)
        self.assertAlmostEqual(r['states']['p']['mean_count'],1-survival,13)
        self.assertAlmostEqual(r['rewards']['photon_count'],.25*(1-survival),13)
        self.assertAlmostEqual(r['rewards']['photon_energy_eV'],.5*(1-survival),13)
        self.assertAlmostEqual(r['checks']['expected_baryons'],2.,13)
        self.assertAlmostEqual(r['checks']['expected_charge'],1.,13)
        self.assertAlmostEqual(sum(v['mean_energy_eV'] for v in r['states'].values()),100.,11)

    def test_zero_length_and_stable_state(self):
        self.assertEqual(solve(fixture(),0.)['initial_state_yield'],1.)
        self.assertEqual(solve(fixture(0.),100.)['rewards']['photon_count'],0.)

    def test_duplicate_fragments_are_counted(self):
        d=fixture();d['initial']='He4';d['states']=[
            dict(id='He4',a=4,z=2,energy_eV=100.,rate_per_Mpc=3.,channels=[dict(probability=1.,products={'D':2})]),
            dict(id='D',a=2,z=1,energy_eV=50.,rate_per_Mpc=0.,channels=[])]
        r=solve(d,.5);self.assertAlmostEqual(r['states']['D']['mean_count'],2*(1-math.exp(-1.5)),13)

    def test_sequential_chain_equal_and_unequal_rates(self):
        for second in [2.,5.]:
            d=fixture();d['states'][0]['channels'][0]['products']={'middle':1}
            d['states'].append(dict(id='middle',a=2,z=1,energy_eV=100.,rate_per_Mpc=second,channels=[dict(probability=1.,products={'p':1,'n':1})]))
            r=solve(d,.7);expected=2*.7*math.exp(-2*.7) if second==2 else 2/(second-2)*(math.exp(-2*.7)-math.exp(-second*.7))
            self.assertAlmostEqual(r['states']['middle']['mean_count'],expected,12)

    def test_analytic_line_bins(self):
        r=line_moments(2.,.25,[0,2,4,math.inf])
        np.testing.assert_allclose(r,[.25,.5,.125,.125,.125,.375,0,0],rtol=0,atol=1e-15)
        self.assertAlmostEqual(sum(r[2::2]),r[0]);self.assertAlmostEqual(sum(r[3::2]),r[1])

    def test_two_path_solvers_agree(self):
        a,b=[solve(fixture(),3.,m) for m in ['expm','ivp']]
        for key in a['rewards']:self.assertAlmostEqual(a['rewards'][key],b['rewards'][key],11)

    def test_distinct_branches_and_conditional_lines(self):
        d=fixture();d['states'][0]['channels']=[
            dict(probability=.3,products={'p':1,'n':1},photons=[dict(mean_energy_eV=2.,probability=.5)]),
            dict(probability=.7,products={'stable_D':1})]
        d['states'].append(dict(id='stable_D',a=2,z=1,energy_eV=100.,rate_per_Mpc=0.,channels=[]))
        r=solve(d,1.);interactions=1-math.exp(-2.)
        self.assertAlmostEqual(r['states']['stable_D']['mean_count'],.7*interactions,13)
        self.assertAlmostEqual(r['states']['p']['mean_count'],.3*interactions,13)
        self.assertAlmostEqual(r['rewards']['photon_count'],.15*interactions,13)

    def test_malformed_or_unphysical_input_fails(self):
        for change in ['rate','prob','missing','multiplicity','charge','bin','energy']:
            d=fixture()
            if change=='rate':d['states'][0]['rate_per_Mpc']=float('nan')
            elif change=='prob':d['states'][0]['channels'][0]['probability']=.5
            elif change=='missing':d['states'][0]['channels'][0]['products']={'absent':1}
            elif change=='multiplicity':d['states'][0]['channels'][0]['products']['p']=.5
            elif change=='charge':d['states'][0]['channels'][0]['products']={'n':2}
            elif change=='bin':d['photon_bin_edges_eV']=[0,4,2,None]
            else:d['energy_scale_eV']=0
            with self.subTest(change=change),self.assertRaises(ValueError):solve(d,1.)


if __name__=='__main__':unittest.main()
