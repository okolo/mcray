"""PD-only coefficient adapters I (independent), M (native), C (official).

The reduction to one energy per isotope is valid only without decay/losses:
all nuclear products retain the injected energy per nucleon. It must not be
used for PD+decay acceptance. No raw photon yields are reinterpreted.
"""
from collections import Counter,defaultdict
import hashlib
import math
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'crpropa'))
from diagnose_network_photo_rates import PhotonIntegral,interpolate,MPC_CM
from run_network_rate_scan import verify

EMITTED=[(1,0),(1,1),(2,1),(3,1),(3,2),(4,2)]


def rows(path):
    return [s.split() for s in path.read_text().splitlines() if s.strip() and not s.startswith('#')]


def identity(a,z):return f'A{a:02}Z{z:02}'


def products(parent,code):
    if len(code)!=6 or not code.isdigit():raise ValueError('Invalid channel code')
    a,z=parent;result=Counter()
    for n,(da,dz) in zip(map(int,code),EMITTED):
        a-=n*da;z-=n*dz
        if n:result[da,dz]+=n
    if not 1<=a<parent[0] or not 0<=z<=a:raise ValueError('Invalid PD residual')
    result[a,z]+=1
    return (a,z),result


def nearest(gamma):
    # Match C++ round for positive gamma-grid coordinates, not Python ties-to-even.
    index=math.floor((math.log10(gamma)-6)/8*200+.5)
    if not 0<=index<=200:raise ValueError('Photon/branch index outside pinned runtime grid')
    return index


def effective_branches(values):
    # The pinned runtime gives any serialized residual probability to the last
    # branch. Do not renormalize its rounded table into a different control.
    if not values or any(not math.isfinite(v) or v<0 for v in values):raise ValueError('Invalid official branches')
    result=[];remaining=1.
    for i,v in enumerate(values):
        prob=remaining if i==len(values)-1 else min(remaining,v)
        result.append(prob);remaining=max(0.,remaining-prob)
    return result


def ratio_probability(numerator,denominator):
    if not math.isfinite(numerator) or not math.isfinite(denominator) or min(numerator,denominator)<0:raise ValueError('Invalid photon rate')
    return min(1.,numerator/denominator) if denominator else float(numerator>0)


class RawBundle:
    def __init__(self,data,photons):
        for name in ['manifest.json']+(['photon_manifest.json'] if photons else []):verify(data,name)
        self.mass={tuple(map(int,r[:2])):float(r[2]) for r in rows(data/'masses.tsv')}
        self.terminal={tuple(map(int,r)) for r in rows(data/'terminals.tsv')}
        self.full={tuple(map(int,r[:2])):np.array(list(map(float,r[2:]))) for r in rows(data/'total.tsv')}
        self.channels=defaultdict(list)
        for r in rows(data/'channels.tsv'):self.channels[tuple(map(int,r[:2]))].append((r[2],np.array(list(map(float,r[3:])))))
        reached=set();pending=[(56,26)]
        while pending:
            p=pending.pop()
            if p in reached:continue
            reached.add(p)
            if p[0]==1 or p in self.terminal:continue
            if p not in self.full or not self.channels[p] or p not in self.mass:raise ValueError('Missing reachable PD data')
            for c,_ in self.channels[p]:pending.extend(products(p,c)[1])
        self.reached=sorted(reached)
        self.grids={key:np.loadtxt(data/f'grid_{key}.tsv')*1e6 for key in ['talys','external']}
        self.transitions={};self.lines=defaultdict(list);self.exceptions=set()
        if photons:
            self.transitions={tuple(map(int,r[:4])):np.array(list(map(float,r[4:]))) for r in rows(data/'pd_photon_transitions.tsv')}
            for r in rows(data/'pd_photons.tsv'):self.lines[tuple(map(int,r[:4]))].append((float(r[4])*1e6,np.array(list(map(float,r[5:])))))
            self.exceptions={tuple(map(int,r)) for r in rows(data/'pd_photon_exceptions.tsv')}
            for p in self.reached:
                for c,_ in self.channels.get(p,[]):
                    key=p+products(p,c)[0]
                    if key not in self.transitions and key not in self.exceptions:raise ValueError('Missing photon coverage')


class Independent:
    def __init__(self,bundle,reference,energy,order):
        e,n=np.loadtxt(reference/'shared_photons.tsv').T;self.photons=PhotonIntegral(e,n)
        self.bundle,self.energy=bundle,energy;self.nodes={};self.weights={}
        nodes,weights=leggauss(order)
        for label,g in bundle.grids.items():
            lo,hi=np.log(g[:-1]),np.log(g[1:]);self.nodes[label]=np.exp(((hi+lo)[:,None]+(hi-lo)[:,None]*nodes)/2).ravel();self.weights[label]=((hi-lo)[:,None]*weights/2).ravel()
        self.integrands={}
    def fold(self,p,xs,photon=False):
        label='talys' if photon or p[0]>=12 else 'external';key=(p,label)
        if key not in self.integrands:
            energy=self.energy[p] if isinstance(self.energy,dict) else self.energy
            gamma=p[0]*energy/(self.bundle.mass[p]*1e6);x=self.nodes[label]
            self.integrands[key]=self.weights[label]*x*x*1e-27*self.photons(x/(2*gamma))*MPC_CM/(2*gamma*gamma)
        value=float(np.sum(self.integrands[key]*interpolate(self.nodes[label],self.bundle.grids[label],xs)))
        if not math.isfinite(value) or value<0:raise ValueError('Invalid independent convolution')
        return value


def build(data,reference,provider,energy,photons=True,order=32,native=None,edges=None,bundle=None,state_energies=None):
    if not math.isfinite(energy) or energy<=0 or order<2:raise ValueError('Invalid energy/quadrature')
    if provider not in ['I','M','C']:raise ValueError('Unknown coefficient provider')
    b=bundle if bundle is not None else RawBundle(data,photons);verify(reference,'shared_manifest.json')
    injection=energy
    if state_energies is not None:
        if set(state_energies)!=set(b.reached) or any(not math.isfinite(e) or e<=0 for e in state_energies.values()):raise ValueError('Invalid state energies')
    mass=dict(b.mass);native_rates={};native_channels={};native_lines=defaultdict(list)
    totals={};branches=defaultdict(list);c_lines=defaultdict(list)
    if provider=='I':fold=Independent(b,reference,state_energies or energy,order).fold
    if provider=='M':
        if native is None:raise ValueError('M requires a native coefficient export')
        text=native.read_text().splitlines()
        header=next((s for s in text if s.startswith('# E_per_A_eV=')),None)
        if header is None or not any(s.startswith('# complete states=') for s in text):raise ValueError('Incomplete native export')
        meta=dict(v.split('=') for v in header[2:].split())
        if float(meta['E_per_A_eV'])!=energy or meta['photons']!=('On' if photons else 'Off'):raise ValueError('Native settings mismatch')
        for r in rows(native):
            p=tuple(map(int,r[1:3]))
            if r[0]=='S':
                if p in native_rates:raise ValueError('Duplicate native state')
                mass[p]=float(r[3]);native_rates[p]=float(r[4])
                if state_energies is not None and (len(r)!=6 or float(r[5])!=state_energies[p]):raise ValueError('Native state energy mismatch')
            elif r[0]=='C':
                key=p+(r[3],)
                if key in native_channels:raise ValueError('Duplicate native channel')
                native_channels[key]=float(r[4])
            elif r[0]=='P':native_lines[p+tuple(map(int,r[3:5]))].append((float(r[5]),float(r[6])))
            elif r[0] in ['T','D','L'] and state_energies is not None:continue
            else:raise ValueError('Unknown native record')
        if set(native_rates)!=set(b.reached):raise ValueError('Native network differs')
        if set(native_channels)!={p+(c,) for p in b.reached for c,_ in b.channels.get(p,[])}:raise ValueError('Native channel coverage differs')
    if provider=='C':
        for r in rows(reference/'Photodisintegration/rate_Shared.txt'):
            z,n=map(int,r[:2]);totals[n+z,z]=np.array(list(map(float,r[2:])))
        for r in rows(reference/'Photodisintegration/branching_Shared.txt'):
            z,n=map(int,r[:2]);branches[n+z,z].append((r[2].zfill(6),np.array(list(map(float,r[3:])))))
        # Preserve the runtime's serialized masses as well as its tables.
        for r in rows(reference/'nuclear_mass.txt'):
            z,n=map(int,r[:2]);mass[z+n,z]=float(r[2])*2.99792458e8**2/1.602176487e-13
        if photons:
            file=reference/'Photodisintegration/photon_emission_Sha.txt'
            if 'Disabled diagnostic photon emission' in file.read_text():raise ValueError('Official photon table is disabled')
            for r in rows(file):
                z,n,zr,nr=map(int,r[:4]);c_lines[z+n,z,zr+nr,zr].append((float(r[4]),np.array(list(map(float,r[5:])))))
    result=[]
    for p in b.reached:
        energy=state_energies[p] if state_energies is not None else injection
        a,z=p;active=a>1 and p not in b.terminal;channels=[]
        gamma=a*energy/(mass[p]*1e6) if a>1 else None
        if provider=='I':
            rate=fold(p,b.full[p]) if active else 0.;raw=[fold(p,xs) for _,xs in b.channels.get(p,[])]
        elif provider=='M':rate=native_rates[p];raw=[native_channels[p+(c,)] for c,_ in b.channels.get(p,[])]
        else:
            in_range=active and 1e6<gamma<1e14
            if in_range and (p not in totals or p not in branches):raise ValueError('Missing official PD state')
            rate=float(np.interp(math.log10(gamma),np.linspace(6,14,201),totals[p])) if in_range else 0.
            active=in_range
            if active and [c for c,_ in branches[p]]!=[c for c,_ in b.channels[p]]:raise ValueError('Official channel order/coverage differs')
            raw=effective_branches([float(v[nearest(gamma)]) for _,v in branches[p]]) if active else []
        if provider!='C' and active:
            total=math.fsum(raw)
            if (rate>0 and total<=0) or (rate==0 and total>0):raise ValueError('Inconsistent full/exclusive rates')
            probabilities=[r/total for r in raw] if total else [1. if i==len(raw)-1 else 0. for i in range(len(raw))]
        else:probabilities=raw
        if active:
            cache={}
            for (code,_),probability in zip(b.channels[p],probabilities):
                residual,emitted=products(p,code);key=p+residual
                if key not in cache:
                    lines=[]
                    if photons:
                        if provider=='I':
                            denominator=fold(p,b.transitions[key],True) if key in b.transitions else 0.
                            source=[(e,ratio_probability(fold(p,xs,True),denominator)) for e,xs in b.lines[key]]
                        elif provider=='M':
                            source=native_lines[key]
                            if len(source)!=len(b.lines[key]):raise ValueError('Native photon coverage differs')
                        else:
                            residual_gamma=residual[0]*energy/(mass[residual]*1e6)
                            source=[(e,min(1.,float(values[nearest(residual_gamma)]))) for e,values in c_lines[key]]
                            if len(source)!=len(b.lines[key]):raise ValueError('Official photon coverage differs')
                        residual_gamma=residual[0]*energy/(mass[residual]*1e6) if source else 0.
                        lines=[dict(mean_energy_eV=e*residual_gamma,probability=pr) for e,pr in source]
                    cache[key]=lines
                channels.append(dict(code=code,probability=probability,products={identity(*k):v for k,v in emitted.items()},photons=cache[key]))
        result.append(dict(id=identity(a,z),a=a,z=z,energy_eV=a*energy,rate_per_Mpc=rate,channels=channels))
    return dict(schema='mcray-deterministic-transport-v1',initial=identity(56,26),energy_scale_eV=56*injection,
        photon_bin_edges_eV=edges or [],states=result,metadata=dict(scope='PD only; decay/losses/EM Off; fixed E/A exact',provider=provider,
        energy_per_nucleon_eV=injection,photons=photons,order=order if provider=='I' else None,
        probability_policy='CRPropa saturated Bernoulli; zero denominator positive numerator => one',boost='residual gamma; isotropic ultrarelativistic',nuclear_energy_debit=False))
