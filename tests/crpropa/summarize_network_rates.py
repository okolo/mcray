#!/usr/bin/env python3
"""Summarize finite-grid rate discrepancies, including absolute weak-rate errors."""
import argparse
import json
import math
from pathlib import Path


def rates(path):
    result={}
    for line in path.read_text().splitlines():
        if line.startswith('# FAILED'):raise ValueError('Failed rate matrix: '+str(path))
        if not line or line.startswith('#'):continue
        a,z,g,c,r=line.split();key=(int(a),int(z),float(g),c);value=float(r)
        if key in result or not math.isfinite(value) or value<0:raise ValueError('Invalid rate matrix')
        result[key]=value
    return result


def worst(rows,metric):
    values=[r for r in rows if r.get(metric) is not None]
    return max(values,key=lambda r:abs(r[metric])) if values else None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--default',type=Path,required=True)
    p.add_argument('--refined',type=Path)
    p.add_argument('--tight',type=Path)
    p.add_argument('--finest',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();reference=json.loads(a.reference.read_text())
    if reference['failed_cells']:raise ValueError('Reference input contains failed rate cells')
    source=reference['rows'];groups={'full':[],'exclusive':[],'sum':[]};sums={}
    for r in source:
        label='full' if r['channel']=='full' else 'exclusive'
        groups[label].append(r)
        if label=='exclusive':
            key=(r['a'],r['z'],r['gamma'])
            entry=sums.setdefault(key,dict(a=r['a'],z=r['z'],gamma=r['gamma'],channel='sum',independent=0.,official_direct=0.,official_gamma_interpolated=0.,refinement_absolute=0.))
            for field in ['independent','official_direct','official_gamma_interpolated','refinement_absolute']:entry[field]+=r[field]
    # Full-only terminal states have an explicitly zero exclusive sum.
    for r in groups['full']:
        key=(r['a'],r['z'],r['gamma'])
        if key not in sums:sums[key]=dict(a=r['a'],z=r['z'],gamma=r['gamma'],channel='sum',independent=0.,official_direct=0.,official_gamma_interpolated=0.,refinement_absolute=0.)
    groups['sum']=list(sums.values())
    summaries={};matrices={}
    for name,path in [('default',a.default),('refined',a.refined),('tight',a.tight),('finest',a.finest)]:
        if path is None:continue
        matrix=rates(path);matrices[name]=matrix;result={};seen=set()
        for label,rows in groups.items():
            compared=[]
            for r in rows:
                key=(r['a'],r['z'],r['gamma'],r['channel']);seen.add(key)
                value=matrix[key];ref=r.get('independent',0.)
                if not math.isfinite(ref) or ref<0:raise ValueError('Invalid independent reference')
                row=dict(a=key[0],z=key[1],gamma=key[2],channel=key[3],rate=value,independent=ref,
                    absolute=value-ref,relative=value/ref-1 if ref else None,refinement_absolute=r.get('refinement_absolute',0.))
                for field in ['official_direct','official_gamma_interpolated']:
                    v=r.get(field,0.)
                    if not math.isfinite(v):raise ValueError('Nonfinite official control')
                    row[field]=v;row[field+'_relative']=v/ref-1 if ref else None
                    row[field+'_absolute']=v-ref
                row['gamma_interpolation_increment']=(r.get('official_gamma_interpolated',0.)-r.get('official_direct',0.))/ref if ref else None
                compared.append(row)
            resolved=[r for r in compared if r['independent']>.1];weak=[r for r in compared if r['independent']<=.1]
            result[label]=dict(rows=len(compared),resolved=len(resolved),weak=len(weak),
                worst_resolved_relative=worst(resolved,'relative'),worst_weak_absolute=worst(weak,'absolute'),
                worst_weak_positive_relative=worst(weak,'relative'),worst_all_absolute=worst(compared,'absolute'),
                zero_reference_nonzero_mcray=sum(r['independent']==0 and r['rate']!=0 for r in compared),
                nonzero_reference_zero_mcray=sum(r['independent']>0 and r['rate']==0 for r in compared),
                largest_reference_with_zero_mcray=max((r['independent'] for r in compared if r['rate']==0),default=0.),
                worst_official_direct=worst(resolved,'official_direct_relative'),
                worst_official_gamma_interpolated=worst(resolved,'official_gamma_interpolated_relative'),
                worst_gamma_interpolation_increment=worst(resolved,'gamma_interpolation_increment'),
                worst_weak_official_direct_absolute=worst(weak,'official_direct_absolute'),
                worst_weak_official_gamma_interpolated_absolute=worst(weak,'official_gamma_interpolated_absolute'),
                worst_reference_refinement=worst(compared,'refinement_absolute'),
                max_resolved_reference_refinement_relative=max((r['refinement_absolute']/r['independent'] for r in resolved),default=0.))
        if seen!=matrix.keys():raise ValueError('Coverage differs between reference and matrix')
        summaries[name]=result
    convergence={}
    for left,right in [('default','refined'),('refined','tight'),('tight','finest')]:
        if left not in matrices or right not in matrices:continue
        x,y=matrices[left],matrices[right]
        if x.keys()!=y.keys():raise ValueError('Refinement matrices differ in coverage')
        by_group={}
        for label,rows in groups.items():
            compared=[]
            for r in rows:
                k=(r['a'],r['z'],r['gamma'],r['channel']);ref=r.get('independent',0.)
                compared.append(dict(a=k[0],z=k[1],gamma=k[2],channel=k[3],absolute=x[k]-y[k],relative=x[k]/y[k]-1 if y[k] else None,resolved=ref>.1))
            by_group[label]=dict(worst_resolved_relative=worst([r for r in compared if r['resolved']],'relative'),worst_weak_absolute=worst([r for r in compared if not r['resolved']],'absolute'))
        convergence[left+'-to-'+right]=by_group
    result=dict(status='finite-grid numerical diagnostic; no physical acceptance verdict',resolution_threshold_per_Mpc=.1,
        gamma_points=len({r['gamma'] for r in source}),gamma_range=[min(r['gamma'] for r in source),max(r['gamma'] for r in source)],
        settings=summaries,convergence=convergence)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({name:{label:dict(resolved=v['resolved'],weak=v['weak'],max_resolved_relative=v['worst_resolved_relative']['relative'] if v['worst_resolved_relative'] else None,max_weak_absolute_per_Mpc=v['worst_weak_absolute']['absolute'] if v['worst_weak_absolute'] else None) for label,v in groups.items()} for name,groups in summaries.items()},indent=2))


if __name__=='__main__':main()
