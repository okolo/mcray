#!/usr/bin/env python3
"""Bounded deterministic comparison of actual historical transport functions.

This reports numerical errors and convergence, not statistical/physical acceptance.
The historical checkout is a hash-pinned, read-only input. No Monte Carlo ensemble
or external interaction data is needed. All generated files go under --output.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import shlex
import subprocess
import sys

import run as kernels

HERE, ROOT = kernels.HERE, kernels.ROOT


def extract(fork, output):
    expected = json.loads((HERE / 'transport_reference.json').read_text())
    actual = {name: kernels.digest(fork / name) for name in expected}
    if actual != expected:
        raise ValueError('Transport reference differs from the pinned snapshot')
    old_time = kernels.block((fork / 'src/MCPropag/Cosmology.h').read_text(), 'class CosmoTime', True)
    new_time = kernels.block((ROOT / 'src/lib/Cosmology.h').read_text(), 'class CosmoTime', True)
    if kernels.normalized_units(old_time) != kernels.normalized_units(new_time):
        raise ValueError('Shared CosmoTime class is not justified by source identity')
    # Reuse the existing verified SavedMF extraction and original RNG files.
    kernels.extract(fork, output)
    header = (fork / 'src/MCPropag/Particle.h').read_text()
    source = (fork / 'src/MCPropag/Particle.cpp').read_text()
    field_header = (fork / 'src/MCPropag/Deflection3D.h').read_text()
    field_source = (fork / 'src/MCPropag/Deflection3D.cpp').read_text()
    fragments = {
        'PARTICLE_ENUM': kernels.block(header, 'enum ParticleType', True),
        'MASS_INLINE': kernels.block(header, 'inline double Mass() const'),
        'CHARGE_INLINE': kernels.block(header, 'inline int ElectricCharge() const'),
        'BETA': kernels.block(header, 'inline double beta() const'),
        'MASSES': kernels.block(source, 'const double Particle::MassesMeV', True),
        'CHARGES': kernels.block(source, 'const int Particle::ElectricCharges', True),
        'MASS_METHOD': kernels.block(source, 'double Particle::Mass(ParticleType'),
        'CHARGE_METHOD': kernels.block(source, 'double Particle::ElectricCharge(ParticleType'),
        'TRACER': kernels.block(source, 'void Particle::updateTracer()'),
        'POSITION_SCALE': kernels.block(field_header, 'virtual double MinVariabilityScale(const long double *x', True),
        'SAVED_CLASS': kernels.block(field_header, 'class SavedMF :', True),
        'SAVED_METHODS': field_source[field_source.index('    SavedMF::SavedMF('):
                                      field_source.index('    void ComboMF::GetValueGauss')],
        'ROTATION': kernels.block(field_source, 'void Deflection3D::GetRotationRate('),
        'RATE': kernels.block(field_source, 'double Deflection3D::Rate('),
        'PROPAGATE': kernels.block(field_source, 'bool Deflection3D::Propagate('),
    }
    template = (HERE / 'transport_reference.cpp.in').read_text()
    for name, value in fragments.items():
        template = template.replace('@' + name + '@', value)
    if '@' in template:
        raise ValueError('Unresolved reference placeholder')
    (output / 'transport_reference.cpp').write_text(template)
    (output / 'fragment-sha256.json').write_text(json.dumps(
        {key: hashlib.sha256(value.encode()).hexdigest()
         for key, value in fragments.items()}, indent=2) + '\n')
    return actual


def inputs(output):
    """Freeze a small case matrix before compilation or trajectory evaluation."""
    grid = output / 'transport.grid'
    with grid.open('w') as stream:
        # Complete, centered, unequal-spacing grid. Each component is independent
        # of its own coordinate: this analytic field and its interpolant have div B=0.
        for i in range(3):
            for j in range(4):
                for k in range(5):
                    x, y, z = -.2 + .2*i, -.3 + .2*j, -.4 + .2*k
                    field = [1e-9*(1 + .3*math.sin(y) + .1*z),
                             1e-9*(.2*math.cos(x) + .1*z*z), 1e-9*(.5 + .2*x*y)]
                    stream.write(' '.join(format(v, '.17g') for v in (x,y,z,*field)) + '\n')
    starts = [('inside', (0,.02,-.03)), ('lower_face', (-.2,.02,-.03)),
              ('upper_face', (.2,.02,-.03)), ('near_outside', (-.25,.02,-.03)),
              ('far_outside', (-1,.02,-.03)), ('corner', (.2,.3,.4))]
    fields = [('free', 0, 0, (0,0,0), starts[0][1]),
              ('axial', 1, 0, (0,0,1e-9), starts[0][1]),
              ('oblique', 1, 0, (2e-10,3e-10,5e-10), starts[0][1]),
              ('gamma10', 2, 0, (0,0,1), starts[0][1])]
    fields += [('grid_'+tag,3,mode,(0,0,0),position)
               for mode in range(4) for tag,position in starts]
    cases = output / 'transport.cases'
    count = 0
    with cases.open('w') as stream:
        for name,kind,mode,field,position in fields:
            for species in range(3):
                for accuracy in (10,20,40):
                    for outer in (1,4,16):
                        # Identical 1 EeV energy, 0.25 Mpc elapsed c*time, fixed
                        # initial direction. gamma10 replaces E and field in the
                        # driver to retain mild curvature at a finite Lorentz factor.
                        row = (name,kind,mode,species,accuracy,outer,1e18,
                               *field,*position,.6,0,.8,.25,.1)
                        stream.write(' '.join(map(str,row)) + '\n')
                        count += 1
    design = dict(cases=count, particles=['electron','positron','proton'],
                  accuracy=[10,20,40], outer_steps=[1,4,16], duration_Mpc=.25,
                  cosmology=dict(hubble_km_s_Mpc=70,omega_lambda=.7,source_z=.01),
                  time_convention='Shared identical CosmoTime class; common initial t(z=.01); no cosmology-kernel comparison.',
                  interactions=[], observers=[], ensembles=False,
                  numerical_acceptance='None: diagnostic errors and convergence only.',
                  sha256={path.name: kernels.digest(path) for path in (grid,cases)})
    (output / 'inputs.json').write_text(json.dumps(design, indent=2) + '\n')
    return cases, grid, design


def summarize(rows):
    groups = {}
    convergence = {}
    for row in rows:
        key = (row['case'],row['mode'],row['particle'],row['accuracy'],row['outer_steps'])
        if key in groups:
            raise ValueError('Duplicate case in trajectory report')
        for implementation in ('fork','modern'):
            for value in row[implementation].values():
                values = value if isinstance(value,list) else [value]
                if any(not math.isfinite(item) for item in values):
                    raise ValueError('Nonfinite trajectory report')
        groups[key] = row
        # Keep field orientation and boundary-start scenarios distinct. Species
        # still share a maximum; the raw trajectory file retains each particle.
        group_key = (row['case'],row['kind'],row['mode'],row['accuracy'],row['outer_steps'])
        group = convergence.setdefault(group_key, dict(count=0,
            max_position_difference_Mpc=0, max_direction_difference=0, max_beta_difference=0,
            fork_zero_scale_calls=0, fork_max_norm_error=0, modern_max_norm_error=0,
            fork_max_time_error_Mpc=0, modern_max_time_error_Mpc=0))
        group['count'] += 1
        for key,field in (('max_position_difference_Mpc','difference_position_Mpc'),
                          ('max_direction_difference','difference_direction'),
                          ('max_beta_difference','difference_beta')):
            group[key] = max(group[key],row[field])
        group['fork_zero_scale_calls'] += row['fork']['zero_scale_calls']
        for implementation in ('fork','modern'):
            values = row[implementation]
            for key,field in (('max_norm_error','direction_norm_error'),
                              ('max_time_error_Mpc','time_error_Mpc'),
                              ('max_analytic_position_error_Mpc','analytic_position_error_Mpc'),
                              ('max_analytic_direction_error','analytic_direction_error')):
                if field in values:
                    target = implementation+'_'+key
                    group[target] = max(group.get(target,0),values[field])
    # Preserve the complete per-accuracy/per-outer-step rows. No post-hoc
    # selection of a best trajectory or a fitted acceptance threshold.
    return dict(cases=len(rows),
                max_position_difference_Mpc=max(r['difference_position_Mpc'] for r in rows),
                max_direction_difference=max(r['difference_direction'] for r in rows),
                fork_zero_scale_calls=sum(r['fork']['zero_scale_calls'] for r in rows),
                convergence=[dict(case=k[0],kind=k[1],mode=k[2],accuracy=k[3],outer_steps=k[4],**v)
                             for k,v in sorted(convergence.items())],
                scope='Deterministic numerical diagnostics; no physical-equivalence verdict.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-root',type=Path,required=True)
    parser.add_argument('--build-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--command-prefix',default='')
    args = parser.parse_args()
    output,build = args.output.resolve(),args.build_dir.resolve()
    fork = args.fork_root.resolve()
    # Reject outputs in the read-only reference, including indirect compiler writes.
    if output == fork or fork in output.parents:
        parser.error('--output must be outside the historical reference')
    output.mkdir(parents=True,exist_ok=True)
    for name in ('report.json','trajectories.jsonl','manifest.json'):
        (output/name).unlink(missing_ok=True)
    reference = extract(fork,output)
    cases,grid,design = inputs(output)
    prefix = shlex.split(args.command_prefix)
    sources = [output/'transport_reference.cpp', HERE/'transport_compare.cpp',
               ROOT/'src/lib/GridMagneticField.cpp', ROOT/'src/lib/Deflection3D.cpp']
    modern = ['Cosmology.h','Cosmology.cpp','Units.cpp','Particle.h','Particle.cpp',
              'Deflection3D.h','Deflection3D.cpp','GridMagneticField.h','GridMagneticField.cpp']
    metadata = dict(reference_sha256=reference,
                    modern_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    modern_working_tree_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
                    modern_sha256={name:kernels.digest(ROOT/'src/lib'/name) for name in modern},
                    adapter_sha256={p.name:kernels.digest(p) for p in HERE.iterdir() if p.is_file()},
                    linked_core_sha256=kernels.digest(build/'mcray/libmcray.a'),
                    environment=dict(platform=platform.platform(),python=sys.version,command_prefix=prefix),
                    reference_adaptations='Verbatim extracted numerical/particle bodies, audited shared CosmoTime/units, driver-owned field lifetime; nuclear branches reject if called.')
    (output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    kernels.compile_probe(build,output,prefix,sources=sources,executable='fork-transport-probe')
    with (output/'execution.log').open('w') as log:
        result = subprocess.run([*prefix,str(output/'fork-transport-probe'),str(cases),str(grid)],
                                text=True,stdout=subprocess.PIPE,stderr=log,timeout=120,check=True)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    if len(rows)!=design['cases']:
        raise ValueError('Incomplete trajectory report')
    report = summarize(rows)
    # Keep original decimal precision of raw C++ coordinates in the JSONL file.
    (output/'trajectories.jsonl').write_text(result.stdout)
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
