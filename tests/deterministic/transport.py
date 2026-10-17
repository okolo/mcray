"""Independent first-moment branching transport and analytic photon rewards.

Coefficient files describe frozen event rates, normalized branch probabilities,
products and Bernoulli lines. No production library or random generator is used.
The solver supplies means, not history variances or an acceptance verdict.
"""
import math
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.integrate import solve_ivp


def positive(value, name, zero=True):
    if isinstance(value,bool) or not math.isfinite(value) or value<0 or (not zero and value==0):
        raise ValueError('Invalid '+name)
    return float(value)


def line_moments(mean_energy, probability, edges):
    """Mean count/energy for a Bernoulli line uniform on [0,2*mean_energy]."""
    mean_energy=positive(mean_energy,'line energy',False)
    width=positive(2*mean_energy,'line upper energy',False)
    probability=positive(probability,'line probability')
    if probability>1:raise ValueError('Line probability exceeds one')
    result=[probability,probability*mean_energy]
    for lower,upper in zip(edges,edges[1:]):
        lo=min(lower,width);hi=min(upper,width)
        fraction=(hi-lo)/width
        result.extend([probability*fraction,probability*fraction*(hi/2+lo/2)])
    return np.array(result)


def bin_edges(values):
    if not values:return []
    edges=[math.inf if v is None else float(v) for v in values]
    if edges[0]!=0 or edges[-1]!=math.inf or any(not math.isfinite(v) for v in edges[:-1]) or any(a>=b for a,b in zip(edges,edges[1:])):
        raise ValueError('Photon bins must strictly partition [0,infinity), ending in null')
    return edges


def operator(model):
    if model.get('schema')!='mcray-deterministic-transport-v1':raise ValueError('Unknown coefficient schema')
    states=model['states'];ids=[s['id'] for s in states]
    if not ids or len(set(ids))!=len(ids):raise ValueError('Missing/duplicate state')
    index={s:i for i,s in enumerate(ids)}
    if model['initial'] not in index:raise ValueError('Missing initial state')
    edges=bin_edges(model.get('photon_bin_edges_eV',[]))
    labels=['photon_count','photon_energy_eV']
    for i in range(max(0,len(edges)-1)):labels.extend([f'photon_bin_{i}_count',f'photon_bin_{i}_energy_eV'])
    energy_scale=positive(model['energy_scale_eV'],'energy scale',False)
    scales=np.array([energy_scale if k.endswith('energy_eV') else 1. for k in labels])
    n=len(states);matrix=lil_matrix((n+len(labels),n+len(labels)),dtype=float)
    for s in states:
        if not isinstance(s['a'],int) or not isinstance(s['z'],int) or s['a']<0:raise ValueError('Invalid state identity')
        positive(s['energy_eV'],'state energy')
    for parent in states:
        i=index[parent['id']];rate=positive(parent['rate_per_Mpc'],'event rate')
        channels=parent['channels']
        if rate and not channels:raise ValueError('Positive rate without channels')
        if channels and not math.isclose(math.fsum(c['probability'] for c in channels),1.,rel_tol=0,abs_tol=1e-12):
            raise ValueError('Branch probabilities do not sum to one')
        matrix[i,i]-=rate
        for channel in channels:
            prob=positive(channel['probability'],'branch probability')
            if prob>1:raise ValueError('Branch probability exceeds one')
            baryons=charge=0
            for product,multiplicity in channel['products'].items():
                if product not in index:raise ValueError('Missing product state: '+product)
                if not isinstance(multiplicity,int) or isinstance(multiplicity,bool) or multiplicity<1:raise ValueError('Invalid multiplicity')
                state=states[index[product]];baryons+=state['a']*multiplicity;charge+=state['z']*multiplicity
                matrix[index[product],i]+=rate*prob*multiplicity
            if (baryons,charge)!=(parent['a'],parent['z']):raise ValueError('Nonconserving channel A/Z')
            for line in channel.get('photons',[]):
                reward=line_moments(line['mean_energy_eV'],line['probability'],edges)/scales
                for k,value in enumerate(reward):matrix[n+k,i]+=rate*prob*value
    return matrix.tocsr(),index,labels,scales


def solve(model,distance,method='expm'):
    distance=positive(distance,'path length')
    matrix,index,labels,scales=operator(model);n=len(index)
    initial=np.zeros(matrix.shape[0]);initial[index[model['initial']]]=1.
    if distance==0:result=initial
    elif method=='expm':result=expm_multiply(matrix*distance,initial)
    elif method=='ivp':
        solution=solve_ivp(lambda x,y:matrix@y,(0,distance),initial,method='DOP853',rtol=2e-12,atol=2e-14)
        if not solution.success:raise ValueError('Independent path integration failed: '+solution.message)
        result=solution.y[:,-1]
    else:raise ValueError('Unknown path solver')
    if not np.isfinite(result).all() or np.min(result)<-1e-10:raise ValueError('Nonfinite/negative transport solution')
    amounts=result[:n];states=model['states']
    return dict(status='deterministic first moments; no acceptance verdict',method=method,distance_Mpc=distance,
        states={s['id']:dict(mean_count=float(amounts[i]),mean_energy_eV=float(amounts[i]*s['energy_eV'])) for i,s in enumerate(states)},
        rewards={key:float(result[n+i]*scales[i]) for i,key in enumerate(labels)},
        checks=dict(expected_baryons=float(sum(v*s['a'] for v,s in zip(amounts,states))),expected_charge=float(sum(v*s['z'] for v,s in zip(amounts,states))),
                    minimum_raw_scaled_component=float(np.min(result))),
        initial_state_yield=float(amounts[index[model['initial']]]))
