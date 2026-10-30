#!/usr/bin/env python3
"""Compile actual hash-pinned fork kernels and compare with the modern source.

The fork is read-only and optional: normal builds and CI do not require it.
All extracted sources, executables, data and reports go into --output.
"""
import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
import random
import re
import shlex
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FILES = ['src/MCPropag/Deflection3D.cpp', 'src/MCPropag/Deflection3D.h',
         'src/CRbeam.cpp', 'src/MCPropag/Jet.cpp', 'src/MCPropag/Jet.h',
         'src/MCPropag/Randomizer.cpp', 'src/MCPropag/Randomizer.h',
         'src/MCPropag/Units.cpp']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block(text, marker, semicolon=False):
    start = text.index(marker)
    opening = text.index('{', start)
    # The pinned fragments have no braces in literals/comments. Assert the
    # exact input hashes first; this is extraction, not a general C++ parser.
    level = 0
    for index in range(opening, len(text)):
        level += (text[index] == '{') - (text[index] == '}')
        if level == 0:
            return text[start:index+1+(1 if semicolon else 0)]
    raise ValueError('Unterminated reference fragment')


def normalized_units(text):
    text = re.sub(r'/\*.*?\*/|//[^\n]*', '', text, flags=re.S)
    return re.sub(r'\s+', '', text.replace('MCPropag', 'mcray'))


def extract(fork, output):
    expected = json.loads((HERE / 'reference.json').read_text())
    actual = {name: digest(fork / name) for name in FILES}
    if actual != expected:
        raise ValueError('Fork input differs from the pinned snapshot; review before changing reference.json')
    if normalized_units((fork / 'src/MCPropag/Units.cpp').read_text()) != normalized_units((ROOT / 'src/lib/Units.cpp').read_text()):
        raise ValueError('Unit implementations differ: a shared Mpc scale is not justified')
    header = (fork / FILES[1]).read_text()
    source = (fork / FILES[0]).read_text()
    saved = block(header, 'class SavedMF :', semicolon=True)
    methods = source[source.index('    SavedMF::SavedMF('):source.index('    void ComboMF::GetValueGauss')]
    cone = block((fork / 'src/CRbeam.cpp').read_text(), 'if(fJetAngle>1.e-14)')
    for name in ('Randomizer.cpp', 'Randomizer.h', 'Jet.cpp', 'Jet.h'):
        (output / name).write_bytes((fork / 'src/MCPropag' / name).read_bytes())
    # Randomizer.cpp includes Utils.h but uses none of its declarations.
    (output / 'Utils.h').write_text('// Empty compatibility include for the isolated reference RNG.\n')
    wrapper = (HERE / 'reference_wrapper.cpp.in').read_text()
    wrapper = wrapper.replace('@SAVED_CLASS@', saved).replace('@SAVED_METHODS@', methods).replace('@CONE_BLOCK@', cone)
    (output / 'reference.cpp').write_text(wrapper)
    return actual


def compile_probe(build, output, prefix, sources=None, executable='fork-kernel-probe'):
    recipe = build / 'tests/CMakeFiles/grid_magnetic_field_tests.dir/link.txt'
    tokens = shlex.split(recipe.read_text())
    compiler, flags, libraries = tokens[0], [], []
    skip = False
    for token in tokens[1:]:
        if skip:
            skip = False
            continue
        if token == '-o':
            skip = True
        elif token.endswith('.o') or 'libgtest' in token:
            continue
        elif token.endswith(('.a', '.so', '.dylib')):
            libraries.append(str((build/'tests'/token).resolve()) if not Path(token).is_absolute() else token)
        elif token.startswith('-l') or token.startswith('-L'):
            libraries.append(token)
        else:
            flags.append(token)
    if sources is None:
        sources = [output/'reference.cpp', HERE/'compare.cpp',
                   ROOT/'src/lib/GridMagneticField.cpp', ROOT/'src/lib/SpatialGeometry.cpp']
    command = [compiler, *flags, '-std=c++17', '-I'+str(ROOT/'src/lib'), '-I'+str(HERE),
               '-iquote', str(output), *map(str, sources),
               '-o', str(output/executable), *libraries]
    (output / 'compiler-version.txt').write_text(subprocess.check_output([*prefix, compiler, '--version'], text=True))
    (output / 'compile-command.json').write_text(json.dumps(command, indent=2)+'\n')
    with (output/'compile.log').open('w') as log:
        subprocess.run([*prefix, *command], cwd=build/'tests', stdout=log, stderr=subprocess.STDOUT, check=True)


def field_inputs(output):
    datasets = []
    for name, low, high, sizes in [('origin', [0,0,0],[1,2,3],[3,4,5]),
                                   ('centered', [-1,-2,-1],[1,2,2],[4,3,5])]:
        grid = output / (name+'.grid')
        with grid.open('w') as stream:
            for i in range(sizes[0]):
                for j in range(sizes[1]):
                    for k in range(sizes[2]):
                        x,y,z = [low[a]+(high[a]-low[a])*v/(sizes[a]-1) for a,v in enumerate((i,j,k))]
                        values = [x,y,z, math.sin(.7*x)+.2*y-.3*z, x*y+.1*z*z, .2*x*z-.4*y]
                        stream.write(' '.join(format(v,'.17g') for v in values)+'\n')
        queries = output/(name+'.queries')
        rng = random.Random(918273)
        with queries.open('w') as stream:
            def emit(tag, fraction):
                position = [low[a]+(high[a]-low[a])*fraction[a] for a in range(3)]
                stream.write(tag+' '+' '.join(format(v,'.17g') for v in position)+'\n')
            for _ in range(1000): emit('interior', [rng.uniform(.001,.999) for _ in range(3)])
            for _ in range(1000): emit('volume', [rng.uniform(-.5,1.5) for _ in range(3)])
            for axis in range(3):
                for value in (-2,-1,-.1,0,1,1.1,2,3):
                    point = [.25,.375,.625]; point[axis]=value
                    emit('face' if value in (0,1) else 'outside_or_seam',point)
        datasets.append((name,grid,queries))
    return datasets


def compare(executable, output, prefix):
    field_rows=[]
    datasets = field_inputs(output)
    (output/'inputs.json').write_text(json.dumps(
        {path.name: digest(path) for _, grid, queries in datasets for path in (grid, queries)},
        indent=2)+'\n')
    for name,grid,queries in datasets:
        for mode in range(4):
            text = subprocess.check_output([*prefix,str(executable),'field',str(grid),str(mode),str(queries)],text=True)
            for line in text.splitlines():
                row=json.loads(line); row['grid']=name; row['mode']=mode; field_rows.append(row)
    groups={}
    for row in field_rows:
        key=(row['grid'],row['mode'],row['category'])
        group=groups.setdefault(key,dict(count=0,bitwise_equal=0,max_absolute=0,max_relative=0,scale_fork=row['scale_fork'],scale_modern=row['scale_modern']))
        group['count']+=1; group['bitwise_equal']+=row['equal']
        group['max_absolute']=max(group['max_absolute'],row['absolute'])
        group['max_relative']=max(group['max_relative'],row['relative'])
    summaries=[dict(grid=k[0],mode=k[1],category=k[2],**value) for k,value in groups.items()]
    (output/'field-queries.json').write_text(json.dumps(field_rows,indent=2)+'\n')
    cone_rows=[]
    # Frozen independent ensembles; no equal-seed trajectory claim.
    for angle in (20.,90.,180.,1e-6*180/math.pi,1e-10*180/math.pi):
        for theta,phi in ((0.,0.),(.4,.7)):
            text=subprocess.check_output([*prefix,str(executable),'cone',str(angle),str(theta),str(phi),'100000'],text=True)
            row=json.loads(text); row.update(angle_degrees=angle,theta=theta,phi=phi)
            cone_rows.append(row)
    comparisons=sum(len(row['metrics']) for row in cone_rows)
    cutoff=statistics.NormalDist().inv_cdf(1-.05/(2*comparisons))
    for row in cone_rows:
        for metric in row['metrics']:
            se=math.sqrt(metric['variance_fork']/row['n']+metric['variance_modern']/row['n'])
            difference=metric['mean_modern']-metric['mean_fork']
            metric.update(difference=difference,standard_error=se,
                          z=difference/se if se else None,
                          simultaneous_95_interval=[difference-cutoff*se,difference+cutoff*se],
                          excludes_zero=abs(difference)>cutoff*se)
    report=dict(field=summaries,cone=cone_rows,confidence=dict(level=.95,method='Bonferroni normal intervals',comparisons=comparisons,z_cutoff=cutoff),
                scope='Kernel comparison only; no full-cascade physical acceptance. Non-detection does not prove zero systematic error.')
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Field queries:',len(field_rows),'identical:',sum(row['equal'] for row in field_rows))
    for row in summaries:
        if row['max_absolute']>1e-10: print('Field discrepancy:',row)
    for row in cone_rows:
        significant=[m['name'] for m in row['metrics'] if m['excludes_zero']]
        print('Cone:',row['angle_degrees'],row['theta'], 'simultaneous intervals excluding zero:',significant)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fork-root',type=Path,required=True)
    p.add_argument('--build-dir',type=Path,required=True,help='Compatible prebuilt modern Release tree with grid test link recipe')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--command-prefix', default='', help='Optional argv prefix for compiler/probe execution, e.g. an offline container command; parsed without a shell')
    args=p.parse_args(); output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    # A failed rerun must not leave an earlier successful report in place.
    for name in ('report.json', 'inputs.json', 'field-queries.json'):
        (output/name).unlink(missing_ok=True)
    prefix = shlex.split(args.command_prefix)
    hashes=extract(args.fork_root.resolve(),output)
    build = args.build_dir.resolve()
    metadata=dict(reference_sha256=hashes,
                  environment=dict(orchestrator_python=sys.version,orchestrator_platform=platform.platform(),command_prefix=prefix),
                  modern_working_tree_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
                  adapter_sha256={path.name: digest(path) for path in sorted(HERE.iterdir()) if path.is_file()},
                  linked_core_sha256=digest(build/'mcray/libmcray.a'),
                  modern_sha256={str(name.relative_to(ROOT)):digest(name) for name in [ROOT/'src/lib/GridMagneticField.cpp',ROOT/'src/lib/GridMagneticField.h',ROOT/'src/lib/SpatialGeometry.cpp']},
                  modern_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  reference_changes='Namespace alias and interface/logging shims only; extracted field methods and source cone block are unchanged.',
                  cone=dict(histories=100000,seeds=dict(fork=17311,modern=91827)))
    (output/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    compile_probe(build,output,prefix)
    compare(output/'fork-kernel-probe',output,prefix)
if __name__=='__main__': main()
