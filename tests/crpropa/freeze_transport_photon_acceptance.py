#!/usr/bin/env python3
"""Freeze initial-point acceptance references from an audited coupled matrix."""
import argparse
import json
from pathlib import Path
from run_transport_photon_pilots import sha

BACKGROUNDS = ['cmb', 'cmb-inoue']
ISOTOPES = ['A55Z26', 'A54Z26', 'A55Z25', 'A54Z25', 'A53Z25']
BINS = [0., 1e14, 1e15, 1e16, None]


def observables(result):
    states = result['states']
    bound = [(int(k[1:].split('Z')[0]), int(k.split('Z')[1]), v)
             for k, v in states.items() if int(k[1:].split('Z')[0]) > 1]
    out = dict(survival=result['initial_state_yield'],
               bound_count=sum(v['mean_count'] for a, z, v in bound),
               bound_A_sum=sum(a*v['mean_count'] for a, z, v in bound),
               bound_Z_sum=sum(z*v['mean_count'] for a, z, v in bound))
    out['energy:bound'] = sum(v['mean_energy_eV'] for a, z, v in bound)
    for species, isotope in [('proton', 'A01Z01'), ('neutron', 'A01Z00')]:
        out['yield:'+species] = states[isotope]['mean_count']
        out['energy:'+species] = states[isotope]['mean_energy_eV']
    for isotope in ISOTOPES:
        out['yield:'+isotope] = states[isotope]['mean_count']
    rewards = result['rewards']
    for species in ['electron', 'positron', 'neutrinoE', 'neutrinoAE']:
        for prefix in ['yield:', 'energy:']:
            out[prefix+species] = rewards[prefix+species]
    out['yield:photon'] = rewards['photon_count']
    out['energy:photon'] = rewards['photon_energy_eV']
    for i in range(4):
        for kind in ['count', 'energy_eV']:
            out[f'photon_bin_{i}_{kind}'] = rewards[f'photon_bin_{i}_{kind}']
    return out


def freeze(matrix):
    hashes = {}
    def read(case, name):
        directory = matrix/case
        manifest = json.loads((directory/'run-manifest.json').read_text())
        if not manifest['status'].startswith('complete PD+decay'):
            raise ValueError('Incomplete reference '+case)
        settings = manifest['settings']
        if (settings['energy_per_nucleon'] != 3.9810717055349857e18 or
                settings['distance_mpc'] != 1 or settings['pd_photons'] != 'On' or
                settings['decay_photons'] != 'On'):
            raise ValueError('Reference scope mismatch')
        for path in ['run-manifest.json', name]:
            hashes[case+'/'+path] = sha(directory/path)
        if manifest['outputs_sha256'][name] != sha(directory/name):
            raise ValueError('Changed reference '+case+'/'+name)
        return json.loads((directory/name).read_text())
    targets = {}; checks = []
    for bg in BACKGROUNDS:
        for provider, code in [('I', 'independent'), ('M', 'mcray'), ('C', 'crpropa')]:
            base = f'{bg}-{provider}-n65'
            means = observables(read(base, 'path-control.json'))
            envelope = {key: 0. for key in means}
            for case, filename in [(base, 'result.json'),
                                   (f'{bg}-{provider}-n33', 'path-control.json')]:
                other = observables(read(case, filename))
                for key in means:
                    envelope[key] += abs(means[key]-other[key])
            if provider == 'I':
                other = observables(read(f'{bg}-I-quadrature', 'path-control.json'))
                for key in means:
                    envelope[key] += abs(means[key]-other[key])
                means = other
            targets[bg+'/'+code] = dict(means=means, numerical_envelope=envelope)
            if provider == 'M':
                for j in [0, 16, 32, 48, 64]:
                    entry = dict(background=bg, case=base, node=j)
                    for label in ['grid', 'native']:
                        name = f'{label}-{j:03}.tsv'
                        manifest = json.loads((matrix/base/'run-manifest.json').read_text())
                        if sha(matrix/base/name) != manifest['outputs_sha256'][name]:
                            raise ValueError('Changed native coefficient control')
                        entry[label+'_sha256'] = sha(matrix/base/name)
                    checks.append(entry)
    input_hashes = {}
    reference_hashes = {}
    for bg in BACKGROUNDS:
        m = json.loads((matrix/f'{bg}-M-n65/run-manifest.json').read_text())
        for path, digest in m['input_sha256'].items():
            if Path(path).name == 'shared_manifest.json':
                reference_hashes[bg] = digest
            else:
                input_hashes[Path(path).name] = digest
    return dict(schema='transport-photon-acceptance-v1', frozen_date='2026-10-20',
        scope='initial Fe56 PD+decay point; both photon sources; no full-domain Release verdict',
        reference_source_commit='7fc04ae6ce2e76887287ba85562b721d6b31ddf8',
        reference_run='20261019T205338Z-943099d8',
        reference_archive_sha256='c9820d862ca39ec01fed92b580beb5b4399fc459b403dafea48b2a3bb45c35c3',
        crpropa_commit='73fa865b3304113da2059755a83f5939c1039d68',
        crpropa_probe_sha256='2aeea9808a3e1c18a6e77df130509d2d6dd8d45df74bfb5f20a516998c148479',
        simultaneous_confidence=.99, integral_tolerance=.01, spectral_tolerance=.05,
        history_cap=200000, prior_pilot_histories=60000, other_controls_reservation=4000,
        main_batches=48, refinement_batches=10, histories_per_batch=2000,
        seed_base=2036102000, energy_per_nucleon_eV=3.9810717055349857e18,
        distance_mpc=1., backgrounds=BACKGROUNDS, photon_bin_edges_eV=BINS,
        selected_isotopes=ISOTOPES, targets=targets, reference_files_sha256=hashes,
        data_manifest_sha256=input_hashes, prepared_reference_manifest_sha256=reference_hashes,
        coefficient_controls=checks,
        native_controls=[dict(step=.0025, eps=1e-6), dict(step=.00125, eps=1e-6), dict(step=.000625, eps=1e-6)],
        reference_step_fractions=[.03/9, .03/27, .03/81])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    result = freeze(a.matrix.resolve())
    with a.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(sha(a.output), a.output)
