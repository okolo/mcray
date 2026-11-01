#!/usr/bin/env python3
"""Deterministic coupling of historical field transport and sphere stopping.

Reuse the approved models and fixtures; report all differences/convergence.
No source RNG ensemble or physical/statistical acceptance verdict is run.
"""
import argparse
from decimal import Decimal
import json
import math
from pathlib import Path
import shlex
import subprocess

import observer

HERE, ROOT = observer.HERE, observer.ROOT
kernels = observer.kernels


def inputs(output):
    # Reuse the exact existing transport grid generator instead of maintaining
    # a second field prescription. Its other generated cases are not executed.
    fixture=output/'transport-fixture'; fixture.mkdir()
    _,grid,transport_design=observer.transport.inputs(fixture)
    profiles=[('pencil_z',0,0),('cone_z',0,0),('cone_rotated',.4,.7)]
    rays=[]
    critical=math.asin(.01/.1)
    near_rim=math.acos(1-.99*(1-math.cos(math.radians(20))))
    for name,theta,phi in profiles:
        axis=[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)]
        center=[.1*x for x in axis]
        offsets=[('axis',0,0)]
        if name!='pencil_z':
            offsets += [(kind,angle,j*math.pi/2) for kind,angle in
                        [('interior',critical/2),('tangent',critical),('near_rim',near_rim)]
                        for j in range(4)]
        for index,(kind,angle,azimuth) in enumerate(offsets):
            # Fixed geometrical input states, not a replacement source sampler.
            # Use the fork's two rotations to prepare the common direction.
            px=math.sin(angle)*math.cos(azimuth); py=math.sin(angle)*math.sin(azimuth)
            pz=math.cos(angle); y=math.cos(theta)*py-math.sin(theta)*pz
            direction=[math.cos(phi+math.pi/2)*px-math.sin(phi+math.pi/2)*y,
                       math.sin(phi+math.pi/2)*px+math.cos(phi+math.pi/2)*y,
                       math.sin(theta)*py+math.cos(theta)*pz]
            norm=math.sqrt(sum(x*x for x in direction)); direction=[x/norm for x in direction]
            rays.append((name,kind+'_'+str(index),center,direction))
    cases=output/'coupled.cases'
    rows=[(name,ray,field,species,accuracy,*center,*direction)
          for name,ray,center,direction in rays for field in (0,1,2)
          for species in (0,1,2) for accuracy in (10,20,40)]
    cases.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows))
    design=dict(cases=len(rows),initial_directions=len(rays),species=['electron','positron','proton'],
        energy_eV=1e18,accuracy=[10,20,40],fields=['none','oblique constant','existing saved grid'],
        constant_field_G=[2e-10,3e-10,5e-10],constant_scale_Mpc=.1,
        observer=dict(mode='strict inside',radius_Mpc=.01,distance_Mpc=.1),
        source_position_Mpc=[0,0,0],source_profiles=[name for name,_,_ in profiles],
        source_sampling=False,interactions=[],duration_Mpc=.2,rate_path='constant',
        cosmology=dict(hubble_km_s_Mpc=70,omega_lambda=.7,source_z=.01),
        acceptance='None: all errors, hit differences and refinement changes retained.',
        tangent_scope='Floating-point/grazing controls, separate from interior/near-rim rows.',
        reused_grid_sha256=transport_design['sha256'][grid.name],case_sha256=kernels.digest(cases))
    (output/'inputs.json').write_text(json.dumps(design,indent=2)+'\n')
    return cases,grid,design


def extract(fork,output):
    reference=observer.extract(fork,output)
    source=(fork/'src/MCPropag/PropagationEngine.cpp').read_text()
    move=kernels.block(source,'void PropagationEngine::move(')
    wrapper=(HERE/'coupled_reference.cpp.in').read_text().replace('@MOVE@',move)
    if '@' in wrapper: raise ValueError('Unresolved coupled fragment')
    (output/'coupled_reference.cpp').write_text(wrapper)
    return reference


def summarize(rows):
    keys=set(); groups={}; trajectories={}
    for row in rows:
        key=tuple(row[name] for name in ('profile','ray','field','particle','accuracy'))
        if key in keys: raise ValueError('Duplicate coupled case')
        keys.add(key)
        if row['same_hit'] != (row['fork']['hit']==row['modern']['hit']):
            raise ValueError('Inconsistent coupled hit comparison')
        for implementation in ('fork','modern'):
            state=row[implementation]
            for value in state.values():
                values=value if isinstance(value,list) else [value]
                if any(not math.isfinite(x) for x in values): raise ValueError('Nonfinite coupled state')
        groupkey=(row['profile'],row['ray'].rsplit('_',1)[0],row['field'],row['particle'],row['accuracy'])
        group=groups.setdefault(groupkey,dict(cases=0,hit_differences=0,max_position_difference_Mpc=0,
            max_direction_difference=0,max_time_difference_Mpc=0,fork_past_end=0,modern_past_end=0,
            fork_zero_scale_calls=0))
        group['cases']+=1; group['hit_differences']+=not row['same_hit']
        for target,source in [('max_position_difference_Mpc','difference_position_Mpc'),
                              ('max_direction_difference','difference_direction'),
                              ('max_time_difference_Mpc','difference_time_Mpc')]:
            if not math.isfinite(row[source]): raise ValueError('Nonfinite difference')
            group[target]=max(group[target],row[source])
        group['fork_past_end']+=row['fork']['past_end']; group['modern_past_end']+=row['modern']['past_end']
        group['fork_zero_scale_calls']+=row['fork']['zero_scale_calls']
        trajectories.setdefault(key[:-1],{})[row['accuracy']]=row
    refinements=[]
    for key,sequence in sorted(trajectories.items()):
        if set(sequence)!={10,20,40}: raise ValueError('Incomplete coupled refinement sequence')
        for low,high in ((10,20),(20,40)):
            entry=dict(profile=key[0],ray=key[1],field=key[2],particle=key[3],accuracy_pair=[low,high])
            for implementation in ('fork','modern'):
                a,b=sequence[low][implementation],sequence[high][implementation]
                entry[implementation]=dict(hit_changed=a['hit']!=b['hit'],
                    position_change_Mpc=math.sqrt(sum((x-y)**2 for x,y in zip(a['position_Mpc'],b['position_Mpc']))),
                    direction_change=math.sqrt(sum((x-y)**2 for x,y in zip(a['direction'],b['direction']))),
                    time_change_Mpc=abs(a['flight_Mpc']-b['flight_Mpc']))
            refinements.append(entry)
    return dict(cases=len(rows),hit_differences=[r for r in rows if not r['same_hit']],
        max_position_difference_Mpc=max(r['difference_position_Mpc'] for r in rows),
        max_direction_difference=max(r['difference_direction'] for r in rows),
        max_time_difference_Mpc=max(r['difference_time_Mpc'] for r in rows),
        groups=[dict(profile=k[0],ray_kind=k[1],field=k[2],particle=k[3],accuracy=k[4],**v)
                for k,v in sorted(groups.items())],refinements=refinements,
        scope='Fixed-state field/observer/engine coupling diagnostic; no sampled-source or physical-equivalence verdict.')


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
    reference=extract(fork,output); cases,grid,design=inputs(output); prefix=shlex.split(args.command_prefix)
    metadata=dict(reference_sha256=reference,
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        adapter_sha256={p.name:kernels.digest(p) for p in HERE.iterdir() if p.is_file()},
        production_sha256={name:kernels.digest(ROOT/'src/lib'/name) for name in
            ('PropagationEngine.cpp','Deflection3D.cpp','GridMagneticField.cpp','SpatialGeometry.cpp','Particle.cpp')},
        linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'),command_prefix=prefix)
    (output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    sources=[output/'coupled_reference.cpp',HERE/'coupled_compare.cpp',HERE/'coupled_modern.cpp',*[ROOT/'src/lib'/name for name in
        ('PropagationEngine.cpp','Deflection3D.cpp','GridMagneticField.cpp','SpatialGeometry.cpp','Jet.cpp')]]
    kernels.compile_probe(build,output,prefix,sources,'coupled-probe',timeout=240)
    with (output/'execution.log').open('w') as log,(output/'trajectories.jsonl').open('w') as raw:
        subprocess.run([*prefix,str(output/'coupled-probe'),str(cases),str(grid)],stdout=raw,stderr=log,check=True,timeout=180)
    rows=[json.loads(line,parse_float=Decimal) for line in (output/'trajectories.jsonl').read_text().splitlines()]
    if len(rows)!=design['cases']: raise ValueError('Incomplete coupled report')
    report=summarize(rows)
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=float)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('groups','refinements','hit_differences')},indent=2,default=float))
    print('Hit differences:',len(report['hit_differences']))


if __name__=='__main__': main()
