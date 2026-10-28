#!/usr/bin/env python3
"""Compare energy moments of two isolated pair-only TransportCR runs."""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def settings(directory):
    root = ET.parse(directory / 'case.xsw').getroot()
    return {e.get('title'): e.get('value') for e in root.iter()
            if e.get('value') is not None}


def spectrum(path):
    lines = path.read_text().splitlines()
    columns = {int(i)-1: name for i, name in
               re.findall(r'(\d+)=([A-Za-z0-9_]+)', lines[0])}
    if not {'p', 'e', 'ae'} <= set(columns.values()):
        raise ValueError(f'Missing required spectrum columns: {path}')
    energy = {name: 0. for name in columns.values()}
    count = dict(energy)
    for line in lines[1:]:
        values = [float(value) for value in line.split()]
        if not values:
            continue
        if values[0] <= 0 or max(columns) >= len(values):
            raise ValueError(f'Invalid spectrum row: {path}')
        for index, name in columns.items():
            energy[name] += values[index]
            count[name] += values[index]/values[0]
    return energy, count


def moments(directory, config):
    initial, counts = spectrum(directory / 'results/case/step0/all')
    final, _ = spectrum(directory / 'results/case/fin/all')
    if any(value != 0. for name, value in initial.items() if name not in {'p', 'All'}):
        raise ValueError('Only primary protons may be injected')
    if counts['p'] <= 0 or initial['p'] <= 0:
        raise ValueError('A nonzero proton source is required')
    injected_energy = float(config['DeltaSourceEnergy'])*1e6
    normalization = counts['p']*injected_energy
    loss = 1.-final['p']/normalization
    emitted = (final['e']+final['ae'])/normalization
    if loss <= 0 or emitted <= 0:
        raise ValueError('Nonzero pair losses and emitted energy are required')
    return {'loss': loss, 'emitted': emitted,
            'balance_residual_relative': abs(emitted-loss)/loss}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse', type=Path, required=True)
    parser.add_argument('--fine', type=Path, required=True)
    parser.add_argument('--threshold', type=float, default=.002)
    args = parser.parse_args()
    if not 0 < args.threshold < 1:
        parser.error('The relative threshold must be between zero and one')
    coarse, fine = settings(args.coarse), settings(args.fine)
    # TransportCR writes these derived normalizations back to the XML.
    numerical = {'energyResolutionFactor', 'microStep', 'Emax',
                 'SourceIntegrDensity', 'SourceIntegrEnergyDensity'}
    if {k: v for k, v in coarse.items() if k not in numerical} != {
            k: v for k, v in fine.items() if k not in numerical}:
        raise ValueError('Physical settings or lower energy bounds differ')
    for config in (coarse, fine):
        if config['MomentSource'] != '1' or config['InjectionSpectraType'] != '1':
            raise ValueError('Monoenergetic instantaneous sources are required')
        if config['NO_PPP_PRODUCTS'] != 'false':
            raise ValueError('Pair secondary production must be enabled')
        if config['CalcMode'] != '1' or config['CELredshift'] != 'false' or config['NO_PPP'] != 'false':
            raise ValueError('Pair-only short-distance runs without redshift losses are required')
        if any(v != 'true' for k, v in config.items()
               if k.startswith('NO_') and not k.endswith('_PRODUCTS')
               and k not in {'NO_PPP', 'NO_MWB'}):
            raise ValueError('Other interactions must be disabled')
    left, right = moments(args.coarse, coarse), moments(args.fine, fine)
    changes = {name: abs(left[name]/right[name]-1.) for name in ('loss', 'emitted')}
    passed = all(change <= args.threshold for change in changes.values())
    print(json.dumps({'coarse': str(args.coarse), 'fine': str(args.fine),
                      'coarse_moments': left, 'fine_moments': right,
                      'relative_changes': changes, 'threshold': args.threshold,
                      'successive_refinement_passed': passed}, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
