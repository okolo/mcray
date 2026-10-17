"""Independent fixed-gamma decay kernels; retain the 51-node beta sampler."""
from collections import Counter,defaultdict
import math
import numpy as np
from pd_coefficients import rows,verify,identity


def beta_mean(q,me):
    """Mean of the actual piecewise-uniform inverse-CDF approximation."""
    if not math.isfinite(q) or q<=0 or not math.isfinite(me) or me<=0:raise ValueError('Invalid beta endpoint')
    x=np.linspace(1.,1.+q/me,51)
    weights=x[1:]*np.sqrt(np.maximum(0.,x[1:]**2-1))*(1+q/me-x[1:])**2
    if not np.isfinite(weights).all() or weights.sum()<=0:raise ValueError('Degenerate beta CDF')
    electron=float(weights@((x[:-1]+x[1:])/2)/weights.sum()*me)
    return electron,q+me-electron


def kernel(parent,code,mass,me):
    """Products and lepton rewards per incoming E/A, in production order."""
    if len(code)!=5 or not code.isdigit() or code=='00000':raise ValueError('Invalid decay code')
    a,z=parent;original_mass=mass[parent]*1e6;gamma_per_epa=a/original_mass
    leptons=Counter()
    for digit,dz,charged,neutral in [(0,1,'electron','neutrinoAE'),(1,-1,'positron','neutrinoE')]:
        for _ in range(int(code[digit])):
            daughter=(a,z+dz)
            q=(mass[a,z]-mass[daughter])*1e6-me
            ee,en=beta_mean(q,me)
            leptons['yield:'+charged]+=1;leptons['yield:'+neutral]+=1
            leptons['energy:'+charged]+=gamma_per_epa*ee;leptons['energy:'+neutral]+=gamma_per_epa*en
            z+=dz
    ratio=mass[a,z]*1e6/original_mass;products=Counter();stopped=False
    for digit,(da,dz) in enumerate([(4,2),(1,1),(1,0)],2):
        for _ in range(int(code[digit])):
            if a<=da or z<dz or z-dz>a-da:raise ValueError('Invalid compound decay residual')
            products[da,dz]+=1;a-=da;z-=dz
            if a>1 and z in [0,a]:
                products[1,int(z>0)]+=a;stopped=True;break
        if stopped:break
    if not stopped:products[a,z]+=1
    if any(p[0]>1 and p not in mass for p in products):raise ValueError('Missing decay product mass')
    if sum(p[0]*n for p,n in products.items())!=parent[0]:raise ValueError('Decay baryon violation')
    charge=sum(p[1]*n for p,n in products.items())-leptons['yield:electron']+leptons['yield:positron']
    if charge!=parent[1]:raise ValueError('Decay charge violation')
    return dict(products=products,ratio=ratio,leptons=leptons)


class DecayBundle:
    def __init__(self,data,reference,raw,provider,photons):
        verify(data,'decay_manifest.json')
        if photons:verify(data,'photon_manifest.json')
        self.mass=dict(raw.mass);self.channels=defaultdict(list)
        self.me=511000.;self.mpc_seconds=3.0856775807e24/2.99792458e10
        stable={tuple(map(int,r[:2])) for r in rows(data/'decay_stable.tsv')}
        native=defaultdict(list);lines=defaultdict(list)
        for r in rows(data/'decays.tsv'):native[tuple(map(int,r[:2]))].append((r[2],float(r[3])))
        if photons:
            covered={(int(r[0]),int(r[1]),r[2]) for r in rows(data/'decay_photon_channels.tsv')}
            if not {(a,z,c) for (a,z),v in native.items() for c,_ in v}<=covered:raise ValueError('Missing decay photon coverage')
            for r in rows(data/'decay_photons.tsv'):lines[int(r[0]),int(r[1]),r[2]].append((float(r[3])*1e3,float(r[4])))
        if provider=='C':
            for r in rows(reference/'nuclear_mass.txt'):
                z,n=map(int,r[:2]);self.mass[z+n,z]=float(r[2])*2.99792458e8**2/1.602176487e-13
            self.me=9.10938291e-31*2.99792458e8**2/1.602176487e-19
            self.mpc_seconds=3.085677581491367e24/2.99792458e10
            official=defaultdict(list);official_lines={}
            for r in rows(reference/'nuclear_decay.txt'):
                z,n=map(int,r[:2]);code=r[2].zfill(5);official[z+n,z].append((code,float(r[3])))
                if (len(r)-4)%2:raise ValueError('Invalid official photon pairs')
                official_lines[z+n,z,code]=[(float(r[i])*1e3,float(r[i+1])) for i in range(4,len(r),2)]
            for p in raw.reached:
                if p[0]==1:continue # Retained isolated fixtures disable free-neutron decay.
                if {c for c,_ in official[p]}!={c for c,_ in native[p]}:raise ValueError('Official decay coverage differs')
            native=official;lines=official_lines if photons else defaultdict(list)
        for p in raw.reached:
            if p[0]==1:continue
            if not native[p] and p not in stable:raise ValueError('Missing explicit decay stability')
            for code,tau in native[p]:
                if not math.isfinite(tau) or tau<=0:raise ValueError('Invalid decay lifetime')
                k=kernel(p,code,self.mass,self.me)
                if any(product not in raw.reached for product in k['products']):raise ValueError('Decay extends prepared network')
                channel=dict(code=code,tau=tau,kernel=k,photons=lines.get((*p,code),[]))
                for e,prob in channel['photons']:
                    if not math.isfinite(e) or e<=0 or not math.isfinite(prob) or not 0<=prob<=1:raise ValueError('Invalid decay photon line')
                self.channels[p].append(channel)


def energy_bounds(raw,decay,initial=(56,26)):
    """DAG interval propagation of all physical paths, independent of rates."""
    from pd_coefficients import products
    bounds={p:[math.inf,-math.inf] for p in raw.reached};bounds[initial]=[1.,1.]
    order=sorted(raw.reached,key=lambda p:(-p[0],-decay.mass.get(p,0.)))
    position={p:i for i,p in enumerate(order)}
    for p in order:
        transitions=[(products(p,c)[1],1.) for c,_ in raw.channels.get(p,[])]
        transitions.extend((v['kernel']['products'],v['kernel']['ratio']) for v in decay.channels[p])
        for children,factor in transitions:
            for child in children:
                if position[child]<=position[p]:raise ValueError('Non-DAG nuclear kernel')
                if bounds[p][0]!=math.inf:
                    bounds[child][0]=min(bounds[child][0],bounds[p][0]*factor)
                    bounds[child][1]=max(bounds[child][1],bounds[p][1]*factor)
    if any(not math.isfinite(v[0]) or v[0]<=0 or v[1]<v[0] for v in bounds.values()):raise ValueError('Unreachable/invalid energy bounds')
    return bounds
