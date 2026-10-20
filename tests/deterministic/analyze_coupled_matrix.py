#!/usr/bin/env python3
"""Audit an extracted full PD+decay matrix and separate numerical refinements.

This reads frozen results only. It runs no transport, consumes no histories,
does not define an acceptance family, and issues no physical acceptance verdict.
"""
import argparse
from bisect import bisect_left
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET


def read(p):return json.loads(p.read_text())


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def payload(root,name):
    p=root/name
    if not p.resolve().is_relative_to(root.resolve()):raise ValueError('Unsafe payload path')
    return p


def close(a,b):
    if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,rel_tol=2e-13,abs_tol=1e-40):
        raise ValueError(f'Inconsistent stored moments: {a}, {b}')


def flat(result):
    values={'survival':result['initial_state_yield'],**result['rewards']}
    for key,state in result['states'].items():
        for label,v in state.items():values[key+':'+label]=v
    return values


def group(key):
    if key.startswith('A'):return 'isotope_energy' if 'energy' in key else 'isotope_count'
    if key.startswith('photon_bin_'):return 'photon_bins'
    if key.startswith(('yield:','energy:')):return 'leptons'
    return 'survival' if key=='survival' else 'photon_totals'


def difference(coarse,fine,atol,energy):
    x,y=flat(coarse),flat(fine)
    if set(x)!=set(y):raise ValueError('Observable coverage changed')
    rows=[];maxima={}
    for key,u in x.items():
        v=y[key];floor=atol*(56*energy if 'energy' in key else 1.)
        relative=(v-u)/v if abs(v)>floor else None
        row=dict(observable=key,coarse=u,fine=v,difference=v-u,relative_to_fine=relative)
        rows.append(row);label=group(key)
        if relative is not None and (label not in maxima or abs(relative)>abs(maxima[label]['relative_to_fine'])):maxima[label]=row
    return dict(maxima=maxima,observables=rows,scope='observed differences; coordinate floors are not certified error bounds')


def restriction(coarse,fine):
    """Restrict fine populations to coarse nodes, preserving count and mean E.

    This is a nodal consistency diagnostic, not a comparison in final physical
    spectral bins. Negative roundoff populations are retained.
    """
    x=defaultdict(list);y=defaultdict(list)
    for r in coarse['energy_states']:x[r['isotope']].append(r)
    for r in fine['energy_states']:y[r['isotope']].append(r)
    rows=[]
    for isotope,old in x.items():
        nodes=[r['energy_per_nucleon_eV'] for r in old];deposited=[0.]*len(nodes)
        for r in y[isotope]:
            e=r['energy_per_nucleon_eV'];v=r['mean_count'];j=bisect_left(nodes,e)
            guard=32*math.ulp(max(e,nodes[-1]))
            if e<nodes[0]-guard or e>nodes[-1]+guard:raise ValueError('Grid bounds changed')
            if j==0:deposited[0]+=v
            elif j==len(nodes):deposited[-1]+=v
            else:
                weight=(e-nodes[j-1])/(nodes[j]-nodes[j-1]);deposited[j-1]+=v*(1-weight);deposited[j]+=v*weight
        l1=math.fsum(abs(v-r['mean_count']) for v,r in zip(deposited,old));yield_fine=fine['states'][isotope]['mean_count']
        rows.append(dict(isotope=isotope,l1_count_per_primary=l1,fine_yield=yield_fine,
            conditional_l1=l1/yield_fine if yield_fine>len(y[isotope])*fine['atol'] else None))
    return dict(total_l1_count_per_primary=math.fsum(r['l1_count_per_primary'] for r in rows),
        largest_absolute=max(rows,key=lambda r:r['l1_count_per_primary']),isotopes=rows)


def audit(results,output):
    meta=read(results/'manifest.json');commit=meta['commit']
    if not re.fullmatch('[0-9a-f]{40}',commit):raise ValueError('Invalid frozen revision')
    if meta['status']!='passed' or any(s.get('exit_code')!=0 for s in meta['steps']):raise ValueError('Incomplete runner')
    for name,digest in meta['payload_sha256'].items():
        if sha(payload(results,name))!=digest:raise ValueError('Changed archive payload: '+name)
    root=results/'deterministic-coupled';matrix=read(root/'matrix.json')
    expected=[]
    for bg in ['cmb','cmb-inoue']:
        for provider in ['I','M','C']:
            for nodes in [17,33,65]:expected.append(dict(bg=bg,provider=provider,nodes=nodes,order=32,step=.0025,id=f'{bg}-{provider}-n{nodes}'))
        expected.extend([dict(bg=bg,provider='I',nodes=65,order=64,step=.0025,id=f'{bg}-I-quadrature'),
                         dict(bg=bg,provider='M',nodes=65,order=32,step=.00125,id=f'{bg}-M-background')])
    if matrix['cases']!=expected or matrix['completed']!=[c['id'] for c in expected] or matrix['histories']!=0 or not matrix['status'].startswith('complete'):raise ValueError('Unexpected/incomplete full matrix')
    repo=Path(__file__).resolve().parents[2]
    @lru_cache(None)
    def git_sha(name):
        blob=subprocess.check_output(['git','show',commit+':'+name],cwd=repo)
        return hashlib.sha256(blob).hexdigest()
    if matrix['driver_sha256']!=git_sha('tests/deterministic/run_coupled_matrix.py'):raise ValueError('Changed matrix driver')
    manifests={};runs={};controls={};binaries=set();timings=[];invariants=[];slice_count=0
    for case in expected:
        name=case['id'];directory=root/name;m=read(directory/'run-manifest.json');s=m['settings']
        if not m['status'].startswith('complete PD+decay') or m['source_commit']!=commit or m['source_status'] or m['histories']!=0:raise ValueError('Incomplete/changed source case: '+name)
        for key,val in [('provider',case['provider']),('nodes',case['nodes']),('order',case['order']),('background_step',case['step']),
                        ('energy_per_nucleon',3.9810717055349857e18),('distance_mpc',1.),('pd_photons','On'),('decay_photons','On'),('path_control',True)]:
            if s[key]!=val:raise ValueError('Changed case settings: '+name+'/'+key)
        if m['completed_coefficient_slices']!=case['nodes']:raise ValueError('Missing coefficient slices')
        for name_input,digest in m['input_sha256'].items():
            if meta['input_sha256'].get(name_input)!=digest:raise ValueError('Changed input identity')
        for name_tool,digest in m['tool_sha256'].items():
            if git_sha(name_tool)!=digest:raise ValueError('Changed frozen tool: '+name_tool)
        for name_output,digest in m['outputs_sha256'].items():
            if sha(payload(directory,name_output))!=digest:raise ValueError('Changed per-case output')
        r=read(directory/'result.json');control=read(directory/'path-control.json');grid=read(directory/'energy-grid.json')
        expected_slices={f'pd-slice-{j:03}.json' for j in range(case['nodes'])}
        if {p.name for p in directory.glob('pd-slice-*.json')}!=expected_slices:raise ValueError('Missing/extra PD coefficient slice')
        for j in range(case['nodes']):
            model=read(directory/f'pd-slice-{j:03}.json');slice_count+=1
            if len(model['states'])!=len(grid) or {v['id'] for v in model['states']}!=set(grid):raise ValueError('Coefficient isotope coverage differs')
            for state in model['states']:
                nodes=grid[state['id']]['nodes_eV'];close(state['energy_eV'],state['a']*nodes[min(j,len(nodes)-1)])
                if not math.isfinite(state['rate_per_Mpc']) or state['rate_per_Mpc']<0:raise ValueError('Invalid coefficient rate')
                if state['channels']:close(math.fsum(c['probability'] for c in state['channels']),1.)
                for channel in state['channels']:
                    if not 0<=channel['probability']<=1:raise ValueError('Invalid coefficient branch')
                    a=z=0
                    for isotope,multiplicity in channel['products'].items():
                        if isotope not in grid or not isinstance(multiplicity,int) or multiplicity<1:raise ValueError('Invalid coefficient product')
                        a+=int(isotope[1:3])*multiplicity;z+=int(isotope[4:6])*multiplicity
                    if (a,z)!=(state['a'],state['z']):raise ValueError('Coefficient channel violates A/Z')
                    for line in channel.get('photons',[]):
                        if not 0<=line['probability']<=1 or not math.isfinite(line['mean_energy_eV']) or line['mean_energy_eV']<=0:raise ValueError('Invalid coefficient line')
        if r['method']!='BDF' or control['method']!='Radau' or control['rtol']!=r['rtol']/10 or control['atol']!=r['atol']/10:raise ValueError('Missing independent stiff path check')
        if len(r['energy_states'])!=m['energy_states'] or len(r['energy_states'])!=sum(len(g['nodes_eV']) for g in grid.values()):raise ValueError('Incomplete energy states')
        for solution in [r,control]:
            populations=defaultdict(list);energies=defaultdict(list)
            for state in solution['energy_states']:
                isotope=state['isotope'];a=int(isotope[1:3]);e=state['energy_per_nucleon_eV']
                if grid[isotope]['nodes_eV'][state['node']]!=e:raise ValueError('Wrong energy-state coordinate')
                populations[isotope].append(state['mean_count']);energies[isotope].append(state['mean_count']*a*e)
            for isotope,state in solution['states'].items():
                close(math.fsum(populations[isotope]),state['mean_count']);close(math.fsum(energies[isotope]),state['mean_energy_eV'])
            rewards=solution['rewards'];close(rewards['photon_count'],rewards['pd_photon_count']+rewards['decay_photon_count'])
            close(rewards['photon_energy_eV'],rewards['pd_photon_energy_eV']+rewards['decay_photon_energy_eV'])
            for label in ['count','energy_eV']:close(rewards['photon_'+label],math.fsum(rewards[f'photon_bin_{k}_{label}'] for k in range(4)))
            baryons=math.fsum(int(k[1:3])*v['mean_count'] for k,v in solution['states'].items())
            charge=math.fsum(int(k[4:6])*v['mean_count'] for k,v in solution['states'].items())+rewards['yield:positron']-rewards['yield:electron']
            energy=math.fsum(v['mean_energy_eV'] for v in solution['states'].values())+math.fsum(rewards['energy:'+k] for k in ['electron','positron','neutrinoE','neutrinoAE'])
            close(baryons,solution['checks']['expected_baryons']);close(charge,solution['checks']['expected_charge_with_leptons']);close(energy/(56*s['energy_per_nucleon']),solution['checks']['energy_without_nondebited_photons_over_initial'])
            invariants.append(dict(case=name,method=solution['method'],baryon_error=baryons-56,charge_error=charge-26,relative_energy_error=energy/(56*s['energy_per_nucleon'])-1,minimum_raw_scaled_component=solution['checks']['minimum_raw_scaled_component']))
        if case['provider']=='M':binaries.add(m['native_probe_sha256'])
        manifests[name]=m;runs[name]=r;controls[name]=control;timings.append(dict(case=name,seconds=m['seconds'],work=r['work'],control_work=control['work']))
    if len(binaries)!=1:raise ValueError('Native binary changed across cases')
    comparisons={};path={};restrictions={};control_grids={}
    for case in expected:
        name=case['id'];s=manifests[name]['settings'];path[name]=difference(runs[name],controls[name],s['atol'],s['energy_per_nucleon'])
    for bg in ['cmb','cmb-inoue']:
        for provider in ['I','M','C']:
            for lo,hi in [(17,33),(33,65)]:
                key=f'{bg}-{provider}-grid-{lo}-{hi}';a=f'{bg}-{provider}-n{lo}';b=f'{bg}-{provider}-n{hi}';s=manifests[b]['settings']
                comparisons[key]=difference(runs[a],runs[b],s['atol'],s['energy_per_nucleon'])
                control_grids[key]=difference(controls[a],controls[b],s['atol']/10,s['energy_per_nucleon'])
                restrictions[key]=restriction(controls[a],controls[b])
        for provider,label in [('I','quadrature'),('M','background')]:
            a=f'{bg}-{provider}-n65';b=f'{bg}-{provider}-{label}';s=manifests[b]['settings'];key=b+'-refinement'
            comparisons[key]=difference(runs[a],runs[b],s['atol'],s['energy_per_nucleon'])
    for key,comparison in comparisons.items():
        stored=read(root/(key+'.json'));actual={r['observable']:r for r in comparison['observables']}
        for row in stored['observables']:
            if any(row[k]!=actual[row['observable']][k] for k in ['coarse','fine','difference','relative_to_fine']):raise ValueError('Inconsistent archived refinement: '+key)
    cross={}
    for bg in ['cmb','cmb-inoue']:
        reference=controls[bg+'-I-quadrature']
        for label,name in [('M_default',bg+'-M-n65'),('M_refined',bg+'-M-background'),('C',bg+'-C-n65')]:
            s=manifests[name]['settings'];comparison=difference(reference,controls[name],s['atol']/10,s['energy_per_nucleon'])
            for row in comparison['observables']:
                floor=s['atol']/10*(56*s['energy_per_nucleon'] if 'energy' in row['observable'] else 1.)
                row['relative_to_independent']=(row['fine']-row['coarse'])/row['coarse'] if abs(row['coarse'])>floor else None
            cross[bg+'-'+label]=comparison
    suite=ET.parse(results/'ctest.xml').getroot();python_log=(results/'workload.log').read_text()
    if suite.get('tests')!='129' or suite.get('failures')!='0' or suite.get('skipped')!='0' or not re.search(r'Ran 21 tests .*?\n\nOK',python_log,re.S):raise ValueError('Missing passing prerequisite checks')
    summary=dict(status='audited full numerical matrix at one injection point; no physical acceptance verdict',source_commit=commit,
        payload_files_verified=len(meta['payload_sha256']),coefficient_slices_checked=slice_count,cases=22,ctest=129,python_tests=21,histories_generated=0,native_binary_sha256=list(binaries),
        analysis_tool_sha256=sha(Path(__file__)),dependencies=manifests[expected[0]['id']]['dependencies'],
        timings=timings,invariants=invariants,refinements=comparisons,path_controls=path,radau_grid_refinements=control_grids,
        energy_population_restrictions=restrictions,coefficient_comparisons=cross,
        finest_means={name:dict(survival=r['initial_state_yield'],rewards=r['rewards']) for name,r in controls.items() if manifests[name]['settings']['nodes']==65})
    output.mkdir(parents=True,exist_ok=False);(output/'matrix-analysis.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(summary['status']);print('verified',summary['payload_files_verified'],'payload files;',summary['cases'],'cases')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--results',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    audit(args.results.resolve(),args.output.resolve())
