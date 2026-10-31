#!/usr/bin/env python3
"""Pointwise comparison with the actual historical vMF weight, not a Monte Carlo verdict."""
import argparse
import json
import math
from pathlib import Path
import shlex
import subprocess

import run as kernels

HERE, ROOT = kernels.HERE, kernels.ROOT


def inputs(output):
    rows=[]
    # Freeze ordinary comparisons separately from deliberate numeric-defect
    # controls; never hide the old zero branch or nonfinite original values.
    for theta,phi in ((0,0),(.4,.7)):
        axis=(math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta))
        directions=[axis,tuple(-x for x in axis),(1,0,0),(0,1,0),(0,0,1)]
        for kappa in (0,.5e-14,-.5e-14,1,-1,100,-100):
            for direction in directions: rows.append(('ordinary',theta,phi,kappa,*direction))
    for kappa in (1e-14,-1e-14,1e-12,-1e-12):
        for z in (-1,0,1): rows.append(('cancellation',0,0,kappa,math.sqrt(1-z*z),0,z))
    for kappa,z in ((1e6,.99925),(-1e6,-.99925),(1e308,1),(-1e308,-1),(-1000,-1)):
        rows.append(('range',0,0,kappa,math.sqrt(max(0,1-z*z)),0,z))
    path=output/'vmf.points'
    path.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows))
    return path,len(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root',type=Path,required=True)
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--command-prefix',default='')
    args=parser.parse_args()
    fork,build,output=args.fork_root.resolve(),args.build_dir.resolve(),args.output.resolve()
    if output==fork or fork in output.parents: parser.error('Output must be outside the reference')
    if output.exists() and any(output.iterdir()): parser.error('Preserve earlier attempts; output must be empty')
    output.mkdir(parents=True,exist_ok=True)
    reference=kernels.extract(fork,output)
    points,count=inputs(output)
    prefix=shlex.split(args.command_prefix)
    metadata=dict(reference_sha256=reference,input_sha256=kernels.digest(points),
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        production_sha256={name:kernels.digest(ROOT/'src/lib'/name) for name in ('Jet.cpp','Jet.h','SpatialGeometry.cpp','SpatialGeometry.h')},
        adapter_sha256={p.name:kernels.digest(p) for p in HERE.iterdir() if p.is_file()},
        linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'),prefix=prefix)
    (output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    sources=[output/'reference.cpp',HERE/'vmf_compare.cpp',ROOT/'src/lib/Jet.cpp']
    kernels.compile_probe(build,output,prefix,sources,'fork-vmf-weight-probe',timeout=240)
    with (output/'execution.log').open('w') as log:
        result=subprocess.run([*prefix,str(output/'fork-vmf-weight-probe'),str(points)],stdout=subprocess.PIPE,
                              stderr=log,text=True,check=True,timeout=120)
    rows=[json.loads(line) for line in result.stdout.splitlines()]
    if len(rows)!=count: raise ValueError('Incomplete pointwise report')
    if any(not math.isfinite(row['modern']) for row in rows): raise ValueError('Nonfinite modern weight')
    (output/'weights.jsonl').write_text(result.stdout)
    ordinary=[row for row in rows if row['group']=='ordinary']
    report=dict(points=count,ordinary_points=len(ordinary),
        max_ordinary_relative_difference=max(row['relative_difference'] for row in ordinary),
        numerical_controls=[row for row in rows if row['group']!='ordinary'],
        scope='Original weight kernel point diagnostics; no source/cascade statistical acceptance.')
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
