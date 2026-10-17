"""PD+decay energy-state transport with conservative remapping and finite rates."""
import math
from collections import defaultdict
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from pd_coefficients import identity
from transport import bin_edges,line_moments,positive

LEPTONS=['electron','positron','neutrinoE','neutrinoAE']


def projection(nodes,energy):
    """Positive linear-energy weights preserve number and mean E/A exactly."""
    nodes=np.asarray(nodes);positive(energy,'remapped energy',False)
    if len(nodes)==0 or not np.isfinite(nodes).all() or nodes[0]<=0 or np.any(np.diff(nodes)<=0):raise ValueError('Invalid energy grid')
    guard=32*np.finfo(float).eps*max(energy,nodes[-1])
    if energy<nodes[0]-guard or energy>nodes[-1]+guard:raise ValueError('Energy outside proven isotope bounds')
    if energy<=nodes[0]:return [(0,1.)]
    if energy>=nodes[-1]:return [(len(nodes)-1,1.)]
    hi=int(np.searchsorted(nodes,energy));fraction=(energy-nodes[hi-1])/(nodes[hi]-nodes[hi-1])
    return [(hi-1,1-fraction),(hi,fraction)]


def grids(bounds,energy,nodes):
    if isinstance(nodes,bool) or nodes<2:raise ValueError('Need at least two energy nodes')
    return {p:np.unique(np.linspace(lo*energy,hi*energy,nodes)) for p,(lo,hi) in bounds.items()}


def assemble(slices,grid,decay,energy,edges=None,pd_photons=True,decay_photons=True,native_decay=None):
    """slices[j] contains PD coefficients at each isotope's j-th energy node."""
    if not slices:raise ValueError('Missing coefficient slices')
    edges=bin_edges(edges or []);states=[];index={}
    for p,nodes in sorted(grid.items()):
        for j,e in enumerate(nodes):index[p,j]=len(states);states.append((p,j,float(e)))
    n=len(states);labels=[]
    for source in ['photon','pd_photon','decay_photon']:
        labels.extend([source+'_count',source+'_energy_eV'])
    for k in range(max(0,len(edges)-1)):labels.extend([f'photon_bin_{k}_count',f'photon_bin_{k}_energy_eV'])
    for species in LEPTONS:labels.extend(['yield:'+species,'energy:'+species])
    reward_index={s:i for i,s in enumerate(labels)}
    scale=56*energy;scales=np.array([scale if ('energy' in label) else 1. for label in labels])
    row=[];column=[];values=[]
    def put(i,j,v):
        if v:row.append(i);column.append(j);values.append(v)
    def reward(i,label,value):put(n+reward_index[label],i,value/scales[reward_index[label]])
    def photons(i,rate,lines,source):
        for line in lines:
            v=line_moments(line['mean_energy_eV'],line['probability'],edges)*rate
            for label,value in zip(['photon_count','photon_energy_eV',*[x for k in range(max(0,len(edges)-1)) for x in [f'photon_bin_{k}_count',f'photon_bin_{k}_energy_eV']]],v):reward(i,label,value)
            reward(i,source+'_count',v[0]);reward(i,source+'_energy_eV',v[1])
    def transition(i,parent,epa,rate,products,factor,leptons=None):
        positive(rate,'channel rate');positive(factor,'energy ratio',False)
        if sum(p[0]*m for p,m in products.items())!=parent[0]:raise ValueError('Nonconserving baryon kernel')
        leptons=leptons or {}
        charge=sum(p[1]*m for p,m in products.items())-leptons.get('yield:electron',0)+leptons.get('yield:positron',0)
        if charge!=parent[1]:raise ValueError('Nonconserving charge kernel')
        for p,multiplicity in products.items():
            if not isinstance(multiplicity,int) or multiplicity<=0:raise ValueError('Invalid product multiplicity')
            for k,w in projection(grid[p],epa*factor):put(index[p,k],i,rate*multiplicity*w)
        for label,v in leptons.items():reward(i,label,rate*v*(epa if label.startswith('energy:') else 1.))
    coefficients=[{s['id']:s for s in model['states']} for model in slices]
    for p,j,epa in states:
        i=index[p,j];s=coefficients[j][identity(*p)];rate=positive(s['rate_per_Mpc'],'PD rate')
        if not math.isclose(s['energy_eV'],p[0]*epa,rel_tol=2e-15):raise ValueError('Coefficient energy mismatch')
        channels=s['channels']
        if rate and not channels:raise ValueError('Missing PD branches')
        if channels and not math.isclose(math.fsum(c['probability'] for c in channels),1,abs_tol=1e-12,rel_tol=0):raise ValueError('Unnormalized PD branches')
        put(i,i,-rate)
        for c in channels:
            pr=positive(c['probability'],'PD probability')
            if pr>1:raise ValueError('Invalid PD probability')
            products={tuple(map(int,key[1:].split('Z'))):m for key,m in c['products'].items()}
            transition(i,p,epa,rate*pr,products,1.)
            if pd_photons:photons(i,rate*pr,c.get('photons',[]),'pd_photon')
        if p[0]==1:continue
        gamma=p[0]*epa/(decay.mass[p]*1e6)
        rates=[decay.mpc_seconds/(c['tau']*gamma) for c in decay.channels[p]]
        lines=[c['photons'] for c in decay.channels[p]]
        total=math.fsum(rates)
        if native_decay is not None:
            native=native_decay[j];expected={c['code'] for c in decay.channels[p]}
            if set(native['channels'].get(p,{}))!=expected:raise ValueError('Native decay channels differ')
            rates=[native['channels'][p][c['code']] for c in decay.channels[p]]
            full=native['totals'][p];summed=math.fsum(rates)
            if (full>0)!=(summed>0):raise ValueError('Inconsistent native decay totals')
            rates=[full*r/summed for r in rates] if summed else [];total=full
            lines=[native['lines'].get((*p,c['code']),[]) for c in decay.channels[p]]
        put(i,i,-positive(total,'decay rate'))
        for c,rate,source in zip(decay.channels[p],rates,lines):
            k=c['kernel'];transition(i,p,epa,rate,k['products'],k['ratio'],k['leptons'])
            if decay_photons:
                if len(source)!=len(c['photons']):raise ValueError('Native decay photon coverage differs')
                photons(i,rate,[dict(mean_energy_eV=e*gamma,probability=pr) for e,pr in source],'decay_photon')
    matrix=coo_matrix((values,(row,column)),shape=(n+len(labels),n+len(labels))).tocsc()
    initial=np.zeros(matrix.shape[0]);initial[index[(56,26),0]]=1.
    return dict(matrix=matrix,initial=initial,states=states,labels=labels,scales=scales,energy_scale_eV=scale)


def solve(operator,distance,method='BDF',rtol=1e-10,atol=1e-16):
    positive(distance,'path length');positive(rtol,'rtol',False);positive(atol,'atol',False)
    if method not in ['BDF','Radau']:raise ValueError('Need a stiff path solver')
    matrix=operator['matrix'];initial=operator['initial'];n=len(operator['states'])
    if distance:
        solution=solve_ivp(lambda x,y:matrix@y,(0,distance),initial,method=method,jac=matrix,rtol=rtol,atol=atol)
        if not solution.success:raise ValueError('Stiff integration failed: '+solution.message)
        final=solution.y[:,-1];work=dict(nfev=solution.nfev,nlu=solution.nlu,steps=len(solution.t))
    else:final=initial;work={}
    if not np.isfinite(final).all() or np.min(final)<-1e-10:raise ValueError('Nonfinite/negative transport result')
    populations=defaultdict(lambda:dict(mean_count=0.,mean_energy_eV=0.));baryons=charge=nuclear_energy=0.;energy_states=[]
    for i,(p,j,epa) in enumerate(operator['states']):
        count=float(final[i]);e=count*p[0]*epa;s=populations[identity(*p)];s['mean_count']+=count;s['mean_energy_eV']+=e
        baryons+=count*p[0];charge+=count*p[1];nuclear_energy+=e
        energy_states.append(dict(isotope=identity(*p),node=j,energy_per_nucleon_eV=epa,mean_count=count))
    rewards={label:float(final[n+k]*operator['scales'][k]) for k,label in enumerate(operator['labels'])}
    charge+=rewards['yield:positron']-rewards['yield:electron']
    no_photon_energy=nuclear_energy+sum(rewards['energy:'+s] for s in LEPTONS)
    return dict(status='PD+decay deterministic diagnostic; no acceptance verdict',method=method,rtol=rtol,atol=atol,distance_Mpc=distance,
        states=dict(populations),energy_states=energy_states,rewards=rewards,initial_state_yield=populations['A56Z26']['mean_count'],work=work,
        checks=dict(expected_baryons=baryons,expected_charge_with_leptons=charge,
            energy_without_nondebited_photons_over_initial=no_photon_energy/operator['energy_scale_eV'],minimum_raw_scaled_component=float(np.min(final))))
