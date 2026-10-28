#!/usr/bin/env python3
"""Audit all reachable decay records and quantify pinned generator corrections."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from prepare_nuclear_decay_data import verify
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
verify(args.data)
args.output.mkdir(parents=True, exist_ok=False)
masses = {tuple(map(int,line.split()[:2])) for line in (args.data/'masses.tsv').read_text().splitlines()
          if line and not line.startswith('#')}
def records(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'): continue
        z, n, code, tau, *_ = line.split()
        key = (int(z)+int(n), int(z), code.zfill(5))
        if key[:2] not in masses: continue
        if key in result: raise ValueError('Duplicate reachable decay record')
        result[key] = float(tau)
    return result
source = args.data / 'decay_sources'
original = records(source / 'reference/data/nuclear_decay.txt')
corrected = records(source / 'corrected/data/nuclear_decay.txt')
inventory = json.loads((source / 'corrected/inventory.json').read_text())
active, stable = set(), set()
for line in (args.data / 'decay_stable.tsv').read_text().splitlines():
    if line and not line.startswith('#'): stable.add(tuple(map(int, line.split())))
rows = []
for line in (args.data / 'decays.tsv').read_text().splitlines():
    if not line or line.startswith('#'): continue
    a, z, code, tau = line.split(); a, z, tau = int(a), int(z), float(tau)
    key = (a, z, code); active.add((a,z))
    if key not in corrected or abs(tau / corrected[key]-1) > 5e-6:
        raise ValueError('Missing or inconsistent corrected lifetime')
    if tau <= 1e-90: raise ValueError('Artificial unknown-isotope dripping is forbidden')
    rows.append((a, z, code, tau, original.get(key), corrected[key],
                 corrected[key]/original[key] if key in original else None))
masses = {tuple(map(int,line.split()[:2])) for line in (args.data/'masses.tsv').read_text().splitlines()
          if line and not line.startswith('#')}
if active & stable or active | stable != masses: raise ValueError('Incomplete decay/stable-state inventory')
for a,z in masses:
    channels = inventory[f'{a},{z}']
    if len(channels)>1 and abs(sum(d['branch'] for d in channels)-1)>1e-10:
        raise ValueError('Exclusive fractions not normalized before ionization')
with (args.output/'lifetimes.csv').open('w') as out:
    writer=csv.writer(out);writer.writerow(('A','Z','channel','mcray_tau_s','original_CRPropa_tau_s','corrected_CRPropa_tau_s','corrected_to_original_tau'))
    writer.writerows(rows)
changed = sorted({(a,z) for a,z,c,t,o,r,ratio in rows if o is None or abs(r/o-1)>1e-5})
summary={'bound_states':len(masses),'active_parent_states':len(active),'stable_or_inactive_states':len(stable),
         'active_channels':len(rows),'affected_isotopes':changed,'status':'structural checks passed'}
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
