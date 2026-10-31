#!/usr/bin/env python3
"""Compare pinned observer kernels and interaction-free ballistic stopping.

This preliminary diagnostic does not sample a source or declare Monte Carlo
equivalence. All differences, including existing numerical safeguards, remain
visible. Historical sources and modern output stay in separate directories.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess

import transport

HERE,ROOT=transport.HERE,transport.ROOT
kernels=transport.kernels


def extract(fork,output):
    expected=json.loads((HERE/'observer_reference.json').read_text())
    actual={name:kernels.digest(fork/name) for name in expected}
    if actual!=expected:
        raise ValueError('Observer reference differs from the pinned snapshot')
    transport.extract(fork,output)
    particle=(fork/'src/MCPropag/Particle.cpp').read_text()
    source=(fork/'src/MCPropag/Observer.cpp').read_text()
    header=(fork/'src/MCPropag/Observer.h').read_text()
    stack=(fork/'src/MCPropag/ParticleStack.h').read_text()
    fragments={
        'FREE_PROPAGATE':kernels.block(particle,'void Particle::PropagateFreely('),
        'OBSERVER_CLASS':kernels.block(header,'class Observer',True),
        'OBSERVER_CONSTRUCTOR':kernels.block(source,'Observer::Observer(int type,'),
        'OBSERVER_REACHED':kernels.block(source,'bool Observer::weInTarget('),
        'OBSERVER_STEP':kernels.block(source,'cosmo_time Observer::patchDtIfNeed('),
        'CONTINUE':kernels.block(stack,'bool ContinuePropagation(const Particle& aParticle) const\n\t{'),
        'TIME_LEFT':kernels.block((fork/'src/MCPropag/ParticleStack.cpp').read_text(),
                                  'cosmo_time Result::GetMinTravelTimeLeft('),
    }
    # Add only the declaration to the shell; the original function body follows.
    path=output/'transport_reference.cpp'
    text=path.read_text()
    assert text.count('void updateTracer();')==1
    path.write_text(text.replace('void updateTracer();',
                                 'void updateTracer();\n    void PropagateFreely(cosmo_time);'))
    text=(HERE/'observer_reference.cpp.in').read_text()
    for key,value in fragments.items(): text=text.replace('@'+key+'@',value)
    if '@' in text: raise ValueError('Unresolved observer fragment')
    (output/'observer_reference.cpp').write_text(text)
    return actual


def inputs(output):
    rows=[]
    # All geometry is fixed before execution. Adjacent-representable points
    # and huge-radius arithmetic are numerical controls, not physical samples.
    points=[('center',0,(0,0,0)),('interior',0,(.5,0,0)),
            ('surface',0,(1,0,0)),('exterior',0,(2,0,0)),
            ('surface_below',1,(1,0,0)),('surface_above',2,(1,0,0)),
            ('diagonal',0,(.6,0,.8))]
    for radius in (.01,.1):
        for center in ((0,0,0),(.25,-.125,.5)):
            for mode in (1,2):
                for name,point,start in points:
                    for duration in (10,.01,1e-9):
                        rows.append((name,0,mode,point,radius,*center,*start,0,0,1,duration))
    for mode in (1,2):
        for name,start in [('range_interior',(.5,0,0)),('range_surface',(1,0,0)),
                           ('range_exterior',(2,0,0))]:
            rows.append((name,0,mode,0,1e165,0,0,0,*start,0,0,1,10))
    flights=[('birth_inside',1,(0,0,0),(0,0,1),4),
             ('birth_outside',2,(0,0,2),(0,0,1),4),
             ('enter',1,(0,0,-3),(0,0,1),6),
             ('exit',2,(0,0,0),(0,0,1),4),
             ('miss',1,(2,0,-3),(0,0,1),6),
             ('surface_inward',1,(0,0,-1),(0,0,1),4),
             ('surface_outward',1,(0,0,1),(0,0,1),4),
             ('short_remaining',1,(0,0,-1),(0,0,1),1e-5),
             ('no_observer',-1,(0,0,0),(0,0,1),4)]
    for radius in (.01,.1):
        for center in ((0,0,0),(.25,-.125,.5)):
            for name,mode,start,direction,duration in flights:
                rows.append((name,1,mode,0,radius,*center,*start,*direction,duration))
    cases=output/'observer.cases'
    cases.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows))
    design=dict(cases=len(rows),points=sum(row[1]==0 for row in rows),
                flights=sum(row[1]==1 for row in rows),
                cosmology=dict(hubble_km_s_Mpc=70,omega_lambda=.7,source_z=.01),
                interactions=[],source_sampling=False,fields=False,
                source='Common initialized photon states; no source-sampler equivalence claim.',
                acceptance='None: deterministic diagnostic, all discrepancies retained.',
                sha256=kernels.digest(cases))
    (output/'inputs.json').write_text(json.dumps(design,indent=2)+'\n')
    return cases,design


def summarize(rows):
    points=[row for row in rows if row['kind']==0]
    flights=[row for row in rows if row['kind']==1]
    return dict(cases=len(rows),points=len(points),flights=len(flights),
                membership_differences=[row for row in points if row['fork_hit']!=row['modern_hit']],
                old_steps_exceeding_remaining=[row for row in points if row['fork_exceeds_remaining']],
                modern_steps_exceeding_remaining=sum(row['modern_exceeds_remaining'] for row in points),
                ballistic_hit_differences=[row for row in flights if row['fork']['hit']!=row['modern']['hit']],
                ballistic_end_time_overshoots=[row for row in flights if row['fork']['past_end']],
                modern_end_time_overshoots=[row for row in flights if row['modern']['past_end'] or row['variable_rate']['past_end']],
                modern_rate_path_differences=[row for row in flights if not row['modern_paths_equal']],
                max_ballistic_position_difference_Mpc=max(row['difference_position_Mpc'] for row in flights),
                max_ballistic_flight_difference_Mpc=max(row['difference_flight_Mpc'] for row in flights),
                max_modern_rate_path_time_difference_Mpc=max(row['difference_rate_path_time_Mpc'] for row in flights),
                scope='Observer kernel and ballistic diagnostics only; no source or cascade statistical acceptance.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root',type=Path,required=True)
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--command-prefix',default='')
    args=parser.parse_args()
    fork,output,build=args.fork_root.resolve(),args.output.resolve(),args.build_dir.resolve()
    if output==fork or fork in output.parents: parser.error('Output must be outside the reference')
    output.mkdir(parents=True,exist_ok=True)
    for name in ('report.json','observer.jsonl','manifest.json'): (output/name).unlink(missing_ok=True)
    reference=extract(fork,output)
    cases,design=inputs(output)
    prefix=shlex.split(args.command_prefix)
    sources=[output/'observer_reference.cpp',HERE/'observer_compare.cpp',HERE/'observer_modern.cpp',
             ROOT/'src/lib/GridMagneticField.cpp',ROOT/'src/lib/Deflection3D.cpp',
             ROOT/'src/lib/SpatialGeometry.cpp',ROOT/'src/lib/PropagationEngine.cpp']
    metadata=dict(reference_sha256=reference,
        modern_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        modern_working_tree_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        adapter_sha256={p.name:kernels.digest(p) for p in HERE.iterdir() if p.is_file()},
        production_sha256={name:kernels.digest(ROOT/'src/lib'/name) for name in
            ('SpatialGeometry.cpp','SpatialGeometry.h','Particle.cpp','ParticleStack.cpp','PropagationEngine.cpp')},
        linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'),command_prefix=prefix)
    (output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    kernels.compile_probe(build,output,prefix,sources=sources,executable='fork-observer-probe')
    with (output/'execution.log').open('w') as log:
        result=subprocess.run([*prefix,str(output/'fork-observer-probe'),str(cases)],
            stdout=subprocess.PIPE,stderr=log,text=True,check=True,timeout=120)
    rows=[json.loads(line) for line in result.stdout.splitlines()]
    if len(rows)!=design['cases']: raise ValueError('Incomplete observer report')
    (output/'observer.jsonl').write_text(result.stdout)
    report=summarize(rows)
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:value if not isinstance(value,list) else len(value) for key,value in report.items()},indent=2))


if __name__=='__main__': main()
