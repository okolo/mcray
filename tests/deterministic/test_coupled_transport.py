"""Analytic controls for finite decay rates, beta kinematics and energy remap."""
from collections import defaultdict
import math
from types import SimpleNamespace
import unittest
import numpy as np
from coupled_transport import projection,grids,assemble,solve
from decay_kernels import beta_mean,kernel,energy_bounds
from pd_coefficients import identity


def fixture(rate=2.,pd_rate=0.):
    parent=(56,26);daughter=(56,27);pd_child=(55,26);neutron=(1,0)
    mass={parent:56000.,daughter:55998.,pd_child:55000.}
    channels=defaultdict(list);k=kernel(parent,'10000',mass,511000.)
    channels[parent]=[dict(code='10000',tau=1/rate,kernel=k,photons=[(1e6,.25)])] if rate else []
    decay=SimpleNamespace(mass=mass,channels=channels,mpc_seconds=1.)
    # E/A=1 GeV gives gamma_parent=1 in this analytic fixture.
    grid={parent:np.array([1e9]),daughter:np.array([1e9*k['ratio']]),pd_child:np.array([1e9]),neutron:np.array([1e9])}
    states=[dict(id=identity(*p),a=p[0],z=p[1],energy_eV=p[0]*e[0],rate_per_Mpc=0.,channels=[]) for p,e in grid.items()]
    states[0].update(rate_per_Mpc=pd_rate,channels=[dict(probability=1.,products={identity(*pd_child):1,identity(*neutron):1},
        photons=[dict(mean_energy_eV=4e5,probability=.3)])] if pd_rate else [])
    model=dict(states=states,schema='mcray-deterministic-transport-v1',initial='A56Z26',energy_scale_eV=56e9,photon_bin_edges_eV=[0,1e6,2e6,None])
    operator=assemble([model],grid,decay,1e9,[0,1e6,2e6,None]);operator['legacy_fixture']=model
    return operator


class CoupledTest(unittest.TestCase):
    def test_projection_number_energy_and_no_clipping(self):
        nodes=[1.,2.,4.];v=projection(nodes,2.75)
        self.assertAlmostEqual(sum(w for j,w in v),1.)
        self.assertAlmostEqual(sum(nodes[j]*w for j,w in v),2.75)
        self.assertTrue(all(w>=0 for j,w in v))
        with self.assertRaises(ValueError):projection(nodes,4.001)
        with self.assertRaises(ValueError):projection([1,1],1.)

    def test_beta_51_node_mean_uses_uniform_cdf_intervals(self):
        q=3.7;me=.511;ee,en=beta_mean(q,me)
        x=me+np.arange(51)*q/50
        density=x*np.sqrt(x*x-me*me)*(q+me-x)**2
        # Integrate the inverse-CDF sampler in probability space independently;
        # midpoint nodes are deterministic quadrature, not Monte Carlo histories.
        cdf=np.cumsum(density);u=(np.arange(100000)+.5)*cdf[-1]/100000
        expected=float(np.mean(np.interp(u,cdf,x)))
        self.assertAlmostEqual(ee,expected,8);self.assertAlmostEqual(ee+en,q+me,14)
        self.assertAlmostEqual(beta_mean(q*1e6,me*1e6)[0]/1e6,ee,14)
        with self.assertRaises(ValueError):beta_mean(-1,1)

    def test_beta_compound_order_and_energy_partition(self):
        mass={(11,4):11000.,(11,5):10998.,(7,3):6998.,(4,2):4000.}
        k=kernel((11,4),'10100',mass,511000.)
        self.assertEqual(k['products'],{(4,2):1,(7,3):1})
        self.assertAlmostEqual(k['ratio'],10998/11000)
        self.assertAlmostEqual(11*k['ratio']+sum(v for name,v in k['leptons'].items() if name.startswith('energy:')),11.)

    def test_finite_decay_survival_and_parent_boosted_photons(self):
        d=fixture();r=solve(d,.7,rtol=1e-12,atol=1e-18);survival=math.exp(-1.4)
        self.assertAlmostEqual(r['initial_state_yield'],survival,9)
        self.assertAlmostEqual(r['states']['A56Z27']['mean_count'],1-survival,9)
        self.assertAlmostEqual(r['rewards']['decay_photon_count'],.25*(1-survival),9)
        self.assertAlmostEqual(r['rewards']['decay_photon_energy_eV'],.25e6*(1-survival),2)
        self.assertAlmostEqual(r['checks']['expected_charge_with_leptons'],26.,10)
        self.assertAlmostEqual(r['checks']['energy_without_nondebited_photons_over_initial'],1.,12)
        self.assertAlmostEqual(r['rewards']['photon_count'],r['rewards']['decay_photon_count'],13)

    def test_competing_pd_decay_rates(self):
        r=solve(fixture(2.,3.),.7);lost=1-math.exp(-5*.7)
        self.assertAlmostEqual(r['initial_state_yield'],1-lost,9)
        self.assertAlmostEqual(r['states']['A55Z26']['mean_count'],lost*3/5,9)
        self.assertAlmostEqual(r['states']['A56Z27']['mean_count'],lost*2/5,9)

    def test_prompt_finite_rate_and_stiff_solver_control(self):
        d=fixture(1e24)
        for method in ['BDF','Radau']:
            r=solve(d,1.,method)
            self.assertAlmostEqual(r['states']['A56Z27']['mean_count'],1.,9)
            self.assertAlmostEqual(r['checks']['expected_baryons'],56.,9)
        # A finite-time observation still resolves this same very short lifetime.
        r=solve(d,1e-24)
        self.assertAlmostEqual(r['initial_state_yield'],math.exp(-1.),9)

    def test_pd_only_reduction_matches_independent_matrix_exponential(self):
        from transport import solve as pd_solve
        d=fixture(0.,2.);r=solve(d,1.,rtol=1e-12,atol=1e-18);old=pd_solve(d['legacy_fixture'],1.)
        self.assertAlmostEqual(r['initial_state_yield'],old['initial_state_yield'],10)
        self.assertAlmostEqual(r['rewards']['photon_count'],old['rewards']['photon_count'],10)
        self.assertAlmostEqual(r['rewards']['photon_energy_eV']/56e9,old['rewards']['photon_energy_eV']/56e9,14)

    def test_dag_bounds_include_distinct_energy_paths(self):
        mass={(56,26):56000.,(56,27):55998.,(55,26):55000.,(55,27):54997.}
        raw=SimpleNamespace(reached=[*mass,(1,0)],channels=defaultdict(list))
        raw.channels[56,26]=[('100000',None)];raw.channels[56,27]=[('100000',None)]
        decay=SimpleNamespace(mass=mass,channels=defaultdict(list))
        for p in [(56,26),(55,26)]:decay.channels[p]=[dict(kernel=kernel(p,'10000',mass,511000.))]
        bounds=energy_bounds(raw,decay);lo,hi=bounds[55,27]
        self.assertAlmostEqual(lo,54997/55000);self.assertAlmostEqual(hi,55998/56000)
        g=grids(bounds,1e18,5)
        for p,v in g.items():self.assertGreaterEqual(v[0],bounds[p][0]*1e18)


if __name__=='__main__':unittest.main()
