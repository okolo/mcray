#!/usr/bin/env python3
"""Run a bounded standalone tail-quantile and sampler timing control."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from run_secondary_server import CASES, EDGES, run, save
from invariant_diagnostics import summarize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=600)
    args=parser.parse_args()
    if not 0<args.seconds<=1800:parser.error('Standalone checks must fit within 30 minutes')
    binary=args.binary.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).resolve().parents[2]
    started=time.monotonic();deadline=started+args.seconds
    edges=output/'edges.txt';edges.write_text('\n'.join(map(str,EDGES))+'\n')
    paths=[]
    for process,energy in CASES:
        case=process+'-'+format(energy,'.0e');path=output/(case+'.tsv')
        run([binary,'diagnostic',process,format(energy,'.17g'),'1e12','.001',edges,path],
            output/(case+'.log'),deadline,cwd=source/'bin')
        paths.append(path)
        print(case,round(time.monotonic()-started,3),flush=True)
    report=summarize(paths)
    report['environment']={'platform':platform.platform(),'machine':platform.machine(),
        'python':sys.version,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'elapsed_seconds':time.monotonic()-started}
    save(output/'summary.json',report)
    print(report['status'])
    return 0 if report['status']=='within_control' else 1


if __name__=='__main__':sys.exit(main())
