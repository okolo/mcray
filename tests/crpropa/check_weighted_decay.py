#!/usr/bin/env python3
"""Finite-path Mn54 beta competition using the real production kernels.

Engineering method controls only, not the final 99% statistical family.
"""
from collections import Counter
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
from run_transport_photon_acceptance import ENV, ROOT
from run_transport_photon_pilots import sha

sys.path.insert(0, str(ROOT/'tests/deterministic'))
from decay_kernels import beta_mean


def check(bundle, references, native_build, ref_build, guides, output):
    output.mkdir(parents=True, exist_ok=False)
    data = bundle/'nuclei'
    native_mass = {tuple(map(int, r[:2])): float(r[2])*1e6
                   for line in (data/'masses.tsv').read_text().splitlines()
                   if line and not line.startswith('#') for r in [line.split()]}
    results = []
    count, gamma, distance = 256, 1e6, 1.
    for ci, code in enumerate(['mcray', 'crpropa']):
        mass = dict(native_mass)
        if code == 'mcray':
            lifetimes = {r[2]: float(r[3]) for line in (data/'decays.tsv').read_text().splitlines()
                         if line.startswith('54 25 ') for r in [line.split()]}
            me, seconds = 511000., 3.0856775807e24/2.99792458e10
        else:
            lifetimes = {r[2]: float(r[3]) for line in (references/'cmb/nuclear_decay.txt').read_text().splitlines()
                         if line.startswith('25 29 ') for r in [line.split()]}
            for line in (references/'cmb/nuclear_mass.txt').read_text().splitlines():
                if line and not line.startswith('#'):
                    z, n, value = line.split()[:3]
                    mass[int(z)+int(n), int(z)] = float(value)*299792458.**2/1.602176487e-19
            me, seconds = 9.10938291e-31*299792458.**2/1.602176487e-19, 3.085677581491367e22/299792458.
        if set(lifetimes) != {'10000', '01000'}:
            raise ValueError('Changed Mn54 competing-decay fixture')
        rate = sum(1/tau for tau in lifetimes.values())
        event_probability = -math.expm1(-distance*seconds*rate/gamma)
        expected = {}
        for channel, species, neutral, dz in [('10000', 'electron', 'neutrinoAE', 1), ('01000', 'positron', 'neutrinoE', -1)]:
            probability = event_probability/(lifetimes[channel]*rate)
            q = mass[54, 25]-mass[54, 25+dz]-me
            charged_energy, neutral_energy = beta_mean(q, me)
            for name, energy in [(species, charged_energy), (neutral, neutral_energy)]:
                expected['yield:'+name] = probability
                expected['energy:'+name] = probability*gamma*energy
        for fi, fraction in enumerate([.5, .8]):
            name = f'{code}-f{fraction:g}'; target = output/(name+'.tsv')
            common = ['--guide', str(guides/'minus.tsv'), '--a', '54', '--charge', '25',
                      '--distance-mpc', str(distance), '--primaries', str(count),
                      '--seed', str(3240000000+ci*100+fi), '--fraction', str(fraction),
                      '--copies', '1', '--output', str(target)]
            if code == 'mcray':
                argv = [str(native_build/'weighted_nuclear_scenario'), '--data', str(data),
                        '--photon-table', str(references/'cmb/shared_photons.tsv'),
                        '--energy-per-nucleon', str(gamma*native_mass[54, 25]/54)]+common
                env = ENV
            else:
                argv = [str(ref_build/'crpropa_weighted_probe'), '--gamma', str(gamma)]+common
                env = {**ENV, 'CRPROPA_DATA_PATH': str(references/'cmb'),
                       'LD_LIBRARY_PATH': str(ref_build/'reference'), 'DYLD_LIBRARY_PATH': str(ref_build/'reference')}
            with (output/(name+'.log')).open('w') as log:
                subprocess.run(argv, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
            histories = {h: Counter() for h in range(count)}
            complete = False
            for line in target.read_text().splitlines():
                if line == f'# batch=0 primaries={count}':
                    complete = True
                if not line or line.startswith('#'):
                    continue
                batch, h, species, weight, energy = line.split()
                if int(batch) != 0 or int(h) not in histories:
                    raise ValueError('Invalid fixture genealogy')
                weight, energy = float(weight), float(energy)
                if not math.isfinite(weight) or weight <= 0 or not math.isfinite(energy) or energy < 0:
                    raise ValueError('Invalid fixture reward')
                histories[int(h)]['yield:'+species] += weight
                histories[int(h)]['energy:'+species] += energy
            if not complete:
                raise ValueError('Incomplete fixture control; never rerun automatically')
            measured = {}
            for key, exact in expected.items():
                values = [h[key] for h in histories.values()]
                mean, se = statistics.mean(values), statistics.stdev(values)/math.sqrt(count)
                if se == 0 or abs(mean-exact) > 6*se:
                    raise ValueError('Weighted finite-path decay sanity control failed or unresolved: '+name+' '+key)
                measured[key] = dict(mean=mean, standard_error=se, expected=exact)
            results.append(dict(id=name, command=argv, primary_executions=count,
                                output_sha256=sha(target), observables=measured))
    result = dict(status='finite-path real-kernel sanity controls passed; not final acceptance',
                  primary_executions=4*count, controls=results,
                  note='Six empirical SEs is a fixed engineering check, not a calibrated 99% joint confidence claim.')
    (output/'controls.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
