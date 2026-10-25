#!/usr/bin/env python3
"""Seeded zero-bias identity controls against original production adapters."""
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import time
from run_transport_photon_acceptance import PROTOCOL, command as ordinary_command
from run_transport_photon_pilots import sha
from run_rare_transport_pilots import command as weighted_command


def rows(path):
    return [tuple(line.split()) for line in path.read_text().splitlines()
            if line and not line.startswith('#')]


def check(bundle, references, native_build, ref_build, guides, output):
    protocol = json.loads(PROTOCOL.read_text())
    output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(guides, output/'guides')
    results = []
    for bi, bg in enumerate(protocol['backgrounds']):
        case = dict(id=bg, background=bg, process='talys-decay', pd_photons=True, decay_photons=True,
                    length_mpc=protocol['distance_mpc'], energy_per_nucleon_eV=protocol['energy_per_nucleon_eV'])
        for ci, code in enumerate(['mcray', 'crpropa']):
            for level in range(3):
                job = dict(id=f'{bg}-{code}-l{level}', case=case, code=code, level=level, histories=2,
                           profile='identity', target='minus', fraction=0., copies=1,
                           seed=3230000000+bi*100+ci*10+level)
                argv, env = weighted_command(job, protocol, bundle, references, native_build, ref_build, output, True)
                argv += ['--force-pd-adapter', '1']
                target = output/(job['id']+'.tsv')
                ordinary = output/(job['id']+'-ordinary.tsv')
                original, original_env = ordinary_command(job, protocol, bundle, references,
                    native_build/'nuclear_scenario', ref_build/'crpropa_decay_probe', ordinary)
                if code == 'crpropa':
                    original_env = {**original_env, 'LD_LIBRARY_PATH': str(ref_build/'reference'),
                                    'DYLD_LIBRARY_PATH': str(ref_build/'reference')}
                start = time.monotonic()
                for label, cmd, variables in [('weighted', argv, env), ('ordinary', original, original_env)]:
                    with (output/(job['id']+'-'+label+'.log')).open('w') as log:
                        subprocess.run(cmd, check=True, env=variables, stdout=log, stderr=subprocess.STDOUT)
                # Compare numerical values to ignore harmless integer formatting.
                def numeric(path):
                    return [(b, h, s, float(w), float(e)) for b, h, s, w, e in rows(path)]
                if {r[1] for r in numeric(target)} != {'0', '1'} or numeric(target) != numeric(ordinary):
                    raise ValueError('Zero-bias seeded output differs from original adapter: '+job['id'])
                results.append(dict(job=job, primary_executions=4, seconds=time.monotonic()-start,
                    command=argv, ordinary_command=original,
                    output_sha256=sha(target), ordinary_sha256=sha(ordinary)))
    result = dict(status='all seeded zero-bias identity controls passed; not physical acceptance',
                  primary_executions=sum(r['primary_executions'] for r in results), controls=results)
    (output/'controls.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['bundle', 'references', 'native-build', 'reference-build', 'guides', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    print(check(a.bundle.resolve(), a.references.resolve(), a.native_build.resolve(),
                a.reference_build.resolve(), a.guides.resolve(), a.output.resolve())['status'])
