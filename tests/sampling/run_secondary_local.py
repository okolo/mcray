#!/usr/bin/env python3
"""Fixed-size, bounded local spectrum sanity check; not server acceptance."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_secondary_server import CASES, EDGES, REJECTION, run, save, seed
import secondary_statistics as stats
import numpy
import scipy


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--events',type=int,default=10000)
    parser.add_argument('--seconds',type=float,default=1200)
    args=parser.parse_args()
    if not 0<args.seconds<=1800 or not 1000<=args.events<=100000:
        parser.error('Require <=30 minutes and 1000--100000 fixed events per mode/case')
    repo=Path(__file__).resolve().parents[2]
    git=lambda *items:subprocess.check_output(['git','-C',str(repo),*items],text=True).strip()
    if git('diff','--name-only','HEAD','--','src/lib','tests/sampling','tests/CMakeLists.txt'):
        parser.error('Commit numerical source and test tools before freezing a local run')
    binary=args.binary.resolve();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    logs=root/'logs';logs.mkdir();started=time.monotonic();deadline=started+args.seconds
    status={'scope':'fixed-size local sanity check; full server precision remains separate',
            'commit':git('rev-parse','HEAD'),'baseline_commit':REJECTION,'N_per_mode_case':args.events,
            'maximum_seconds':args.seconds,'status':'running',
            'environment':{'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,
                           'platform':platform.platform()}}
    save(root/'status.json',status)
    try:
        run(['git','-C',repo,'archive','--format=tar','--output='+str(root/'source.tar'),'HEAD'],logs/'source.log',deadline)
        shutil.copy2(binary,root/'em_secondary_probe')
        cache=binary.parent.parent/'CMakeCache.txt'
        if cache.is_file():shutil.copy2(cache,root/'CMakeCache.txt')
        manifest=binary.parent/'sampling_reference/manifest.json'
        if manifest.is_file():shutil.copy2(manifest,root/'build-manifest.json')
        edges=root/'edges.txt';edges.write_text('\n'.join(map(str,EDGES))+'\n')
        references={};counts={};seeds={}
        for process,energy in CASES:
            case=process+'-'+format(energy,'.0e');directory=root/case;directory.mkdir()
            common=[process,format(energy,'.17g'),'1e12','.001',edges]
            paths=[]
            for cells,order in ((128,16),(256,24),(512,32)):
                path=directory/('oracle-'+str(cells)+'.tsv');paths.append(path)
                run([binary,'oracle',*common,cells,order,path],logs/(case+'-oracle.log'),deadline,cwd=repo/'bin')
            references[case]=stats.freeze_reference(paths);counts[case]={};seeds[case]={}
            for mode,algorithm in (('old','rejection'),('cdf','cdf'),('current','current')):
                value=seed(case+'local-bounded-20261020'+mode);seeds[case][mode]=value
                stem=directory/mode
                run([binary,'sample',*common,algorithm,args.events,value,str(stem)+'.events',str(stem)+'.tsv'],
                    logs/(case+'-'+mode+'.log'),deadline,cwd=repo/'bin')
                if Path(str(stem)+'.events').stat().st_size!=8+26*args.events:
                    raise RuntimeError('Incomplete local raw event file')
                counts[case][mode]=stats.read_counts([Path(str(stem)+'.tsv')])
            print(case,round(time.monotonic()-started,3),flush=True)
        labels={'old':'Frozen rejection (1cc28f3)','cdf':'Numerical log-CDF control','current':'Analytic bounded rejection','old_minus_current':'Control minus bounded'}
        data={'references':references,'counts':counts,'planned_n':{case:args.events for case in references},
              'complete':True,'sampler_labels':labels,'controls':['old','cdf']}
        save(root/'analysis-input.json',data);save(root/'seeds.json',seeds)
        run([sys.executable,repo/'tests/sampling/analyze_bounded_secondary.py',root/'analysis-input.json',root/'report.json'],
            logs/'analysis.log',deadline)
        report=json.loads((root/'report.json').read_text())
        status['status']=report['status'];status['main_events']=args.events*3*len(CASES)
        status['binary_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
    except Exception as error:
        status.update(status='execution_error',error=str(error))
        raise
    finally:
        status['elapsed_seconds']=time.monotonic()-started;save(root/'status.json',status)
        manifest={}
        for path in root.rglob('*'):
            if path.is_file():manifest[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
        save(root/'sha256.json',manifest)
    print(json.dumps(status))
    return {'within_target':0,'deviation':2,'inconclusive':3}.get(status['status'],1)


if __name__=='__main__':sys.exit(main())
