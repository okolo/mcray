#!/usr/bin/env python3
"""Prepare an isolated TransportCR XML run without modifying its checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

p = argparse.ArgumentParser()
p.add_argument('--transportcr', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--process', choices=('none', 'psb', 'gzk', 'pair'), required=True)
p.add_argument('--a', type=int, default=56)
p.add_argument('--energy-per-nucleon', type=float, required=True)
p.add_argument('--distance-mpc', type=float, required=True)
p.add_argument('--resolution', type=int, default=4)
p.add_argument('--step-mpc', type=float, default=0.1)
p.add_argument('--grid-span-decades', type=float, default=3.)
p.add_argument('--sophia-tables', type=Path,
               help='Copy an external SOPHIA table set into the isolated run')
a = p.parse_args()
if not 0 <= a.a <= 56 or a.resolution < 1 or min(a.energy_per_nucleon, a.distance_mpc, a.step_mpc, a.grid_span_decades) <= 0:
    p.error('Invalid physical/grid parameters')
a.output.mkdir(parents=True, exist_ok=False)
root = ET.parse(a.transportcr / 'bin/switches.xsw').getroot()
params = {e.get('title'): e for e in root.iter() if e.get('value') is not None}
values = {
    'CalcMode': '1', 'MomentSource': '1', 'CELredshift': 'false',
    'L_short_distance_test': str(a.distance_mpc), 'sdtTauParticle': '-1',
    'sdtDiffusionEmax': '0', 'LeakyBoxTest': 'false', 'shortDistanceTestOutputNo': '4',
    'InjectionSpectraType': '1', 'AutoNucleiRates': 'false',
    'DeltaSourceEnergy': str(a.energy_per_nucleon * max(a.a, 1) / 1e6),
    'SourceNormMode': '0', 'SourceAbsNorm': '1', 'EmaxMode': '2',
    'Emin': str(a.energy_per_nucleon / 10**(a.grid_span_decades/2)),
    'Emax': str(a.energy_per_nucleon * 10**(a.grid_span_decades/2)),
    'energyResolutionFactor': str(a.resolution), 'microStep': str(a.step_mpc),
    'NO_MWB': 'false', 'T0_CMB': '2.7255', 'IRO_model': '0', 'Radio_model': '0',
    'NO_PHOTODISINTEGRATION': str(a.process != 'psb').lower(),
    'NO_PI': str(a.process != 'gzk').lower(),
    'NO_PPP': str(a.process != 'pair').lower(),
    'NO_PI_NUCLEI': 'true', 'COEF_TEST_ON': 'true', 'Zmax': '0.01',
    'H_in_km_s_Mpc': '67.8', 'Lv': '0.692',
    'MaxNucleiA': str(56 if a.a > 1 else 1), 'B0': '0',
    'NO_PI_PRODUCTS': 'false', 'NO_PPP_PRODUCTS': 'false',
}
# All other processes remain disabled, including decay and secondary cascades.
for key in params:
    if key.startswith('NO_') and key not in values:
        values[key] = 'true'
for key in ('eRate', 'aeRate', 'phRate', 'enRate', 'aenRate', 'mnRate', 'amnRate',
            'tnRate', 'atnRate', 'pRate', 'nRate'):
    values[key] = '0'
for mass in range(2, 57):
    values[f'A{mass:02}Rate'] = '0'
values['nRate' if a.a == 0 else 'pRate' if a.a == 1 else f'A{a.a:02}Rate'] = '1'
for key, value in values.items():
    if key not in params:
        raise RuntimeError(f'Missing TransportCR XML parameter: {key}')
    params[key].set('value', value)
for key in params:
    if key in ('Bgauss', 'Bfield'):
        params[key].set('value', '0')
xml = ET.tostring(root, encoding='unicode')
(a.output / 'case.xsw').write_text('<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!DOCTYPE dir SYSTEM "switches.dtd">\n' + xml, encoding='utf-8')
shutil.copyfile(a.transportcr / 'bin/switches.dtd', a.output / 'switches.dtd')
reference_tables = (a.transportcr / 'bin/tables').resolve()
if a.sophia_tables is None:
    (a.output / 'tables').symlink_to(reference_tables, target_is_directory=True)
else:
    source = a.sophia_tables.expanduser().resolve()
    density = str(10 * a.resolution)
    files = [source / 'p', source / 'n']
    files += [source / density / f'{primary}-{secondary}'
              for primary in ('p', 'n')
              for secondary in ('p', 'n', 'e', 'ae', 'ph', 'en', 'aen', 'mn', 'amn')]
    for file in files:
        if not file.is_file():
            raise RuntimeError(f'Missing SOPHIA table: {file}')
        if file.parent.name == density:
            header = file.read_text().splitlines()[0].split()
            if len(header) != 5 or int(header[3]) != int(density):
                raise RuntimeError(f'Incompatible SOPHIA density: {file}')
    tables = a.output / 'tables'
    tables.mkdir()
    for child in reference_tables.iterdir():
        if child.name != 'sophia2':
            (tables / child.name).symlink_to(child, target_is_directory=child.is_dir())
    destination = tables / 'sophia2'
    (destination / density).mkdir(parents=True)
    hashes = {}
    for file in files:
        relative = file.relative_to(source)
        shutil.copyfile(file, destination / relative)
        hashes[str(relative)] = hashlib.sha256((destination / relative).read_bytes()).hexdigest()
    (a.output / 'sophia-manifest.json').write_text(json.dumps(
        {'source': str(source), 'bins_per_decade': int(density), 'sha256': hashes},
        indent=2) + '\n', encoding='utf-8')
