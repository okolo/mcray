#!/usr/bin/env python3
"""Independently recount saved binary events and verify their physical support."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np


def audit(root, seconds=600):
    deadline=time.monotonic()+seconds
    def check_time():
        if time.monotonic()>deadline:raise TimeoutError('Raw audit deadline reached')
    dtype=np.dtype([('s','<f8'),('a','<f8'),('b','<f8'),('type','u1'),('flags','u1')])
    assert dtype.itemsize==26
    result={'events':0,'files':{},'hashes_checked':0}
    for raw in sorted(root.rglob('*.events')):
        check_time();counts=raw.with_suffix('.tsv')
        if not counts.is_file():raise ValueError('Missing completed counts: '+str(raw))
        metadata=dict(item.split('=') for item in counts.read_text().splitlines()[0][2:].split())
        n=int(metadata['N'])
        with raw.open('rb') as source:
            if source.read(8)!=b'MCRAYEM1' or raw.stat().st_size!=8+26*n:raise ValueError('Invalid raw size/magic')
        records=np.memmap(raw,dtype=dtype,mode='r',offset=8,shape=(n,))
        a=records['a'].astype(np.longdouble);b=records['b'].astype(np.longdouble);s=records['s'].astype(np.longdouble)
        if not np.all(np.isfinite(a)&np.isfinite(b)&np.isfinite(s)&(a>0)&(a<1)&(b>0)&(b<1)):
            raise ValueError('Nonphysical raw value')
        error=float(np.max(np.abs(a+b-1)));tolerance=32*np.finfo(float).eps
        pair=metadata['process']=='PP'
        if pair:
            if np.any(s<4) or np.any(~np.isin(records['type'],[0,1])):raise ValueError('Invalid PP types/invariant')
            lower=2/(s*(1+np.sqrt(1-4/s)));upper=np.longdouble(.5)
        else:
            if np.any(s<=1) or np.any(records['type']!=2):raise ValueError('Invalid ICS types/invariant')
            lower=np.longdouble(metadata['Emin_eV'])/np.longdouble(metadata['E_eV']);upper=1-1/s
        violation=float(max(0,np.max(lower-a),np.max(a-upper)))
        if error>tolerance or violation>tolerance:raise ValueError('Energy/support violation')
        expected={line.split()[0]:list(map(int,line.split()[1:])) for line in counts.read_text().splitlines() if not line.startswith('#')}
        oracle=max(raw.parent.glob('oracle-*.tsv'),key=lambda p:int(p.stem.split('-')[-1]))
        electron=records['type']==0
        views={'gamma':a,'low':a,'recoil':b,'high':b,'electron':np.where(electron,a,b),'positron':np.where(electron,b,a)}
        for line in oracle.read_text().splitlines():
            if line.startswith('#'):continue
            check_time();fields=line.split();key,view=fields[:2];lo,hi=map(np.longdouble,fields[3:5])
            def inside(x):return ((x>=lo)&((x<hi)|((hi==1)&(x==1)))).astype(np.int8)
            multiplicity=inside(a)+inside(b) if view=='inclusive' else inside(views[view])
            actual=np.bincount(multiplicity,minlength=3).tolist()
            if actual!=expected.pop(key):raise ValueError('Raw bin mismatch: '+str(raw)+' '+key)
        if expected:raise ValueError('Missing raw observable')
        if metadata['sampler']=='current' and (np.any(records['flags']) or int(metadata.get('bounded_fallbacks',0))):
            raise ValueError('Bounded run recorded a fallback/violation')
        result['events']+=n
        result['files'][str(raw.relative_to(root))]=dict(N=n,max_energy_error=error,max_support_error=violation,
                                                     flagged_events=int(np.count_nonzero(records['flags'])))
    if not result['files']:raise ValueError('No raw event files')
    manifest=root/'sha256.json'
    if manifest.is_file():
        for name,expected in json.loads(manifest.read_text()).items():
            h=hashlib.sha256()
            with (root/name).open('rb') as source:
                for chunk in iter(lambda:source.read(1024*1024),b''):check_time();h.update(chunk)
            if h.hexdigest()!=expected:raise ValueError('Checksum mismatch: '+name)
            result['hashes_checked']+=1
    result['all_raw_counts_reproduced']=True
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--seconds',type=float,default=600)
    args=parser.parse_args()
    if not 0<args.seconds<=1800:parser.error('Require 0 < seconds <= 1800')
    result=audit(args.root,args.seconds)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Verified',result['events'],'events and',result['hashes_checked'],'hashes')
