#!/usr/bin/env python3
"""Inventory one-step GZK spectators in a prepared nuclear data bundle.

This is not a full combined-network closure computation. Missing treatment
is reported, never inferred from the photodisintegration terminal whitelist.
"""
import argparse
import hashlib
import json
from pathlib import Path


def identities(path):
    result = set()
    for number, line in enumerate(path.read_text().splitlines(), 1):
        row = line.split('#', 1)[0].split()
        if not row:
            continue
        if len(row) < 2:
            raise ValueError(f'{path}:{number}: missing A,Z')
        a, z = map(int, row[:2])
        if a < 2 or not 0 < z <= a:
            raise ValueError(f'{path}:{number}: invalid bound A,Z')
        result.add((a, z))
    return result


def audit(data):
    files = ['masses.tsv', 'decays.tsv', 'decay_stable.tsv', 'total.tsv']
    masses, decays, stable, pd = [identities(data / f) for f in files]
    missing = {}
    for a, z in sorted(masses):
        for target, dz, count in [('p', 1, z), ('n', 0, a-z)]:
            if not count:
                continue
            residual = (a-1, z-dz)
            if residual in {(1, 0), (1, 1)}:
                continue
            gaps = [kind for kind, available in
                    [('mass', masses), ('decay', decays | stable), ('PD', pd)]
                    if residual not in available]
            if gaps:
                entry = missing.setdefault(residual, {
                    'A': residual[0], 'Z': residual[1], 'missing': gaps, 'from': []})
                entry['from'].append({'A': a, 'Z': z, 'target': target})
    return {
        'scope': 'one-step spectators from prepared bound-state inventory; not full closure',
        'bundle': str(data.resolve()), 'bound_states': len(masses),
        'input_sha256': {f: hashlib.sha256((data/f).read_bytes()).hexdigest() for f in files},
        'missing_states': list(missing.values())}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = audit(args.data)
    with args.output.open('x') as output:
        json.dump(report, output, indent=2)
        output.write('\n')
    print(f"{report['bound_states']} bound states; "
          f"{len(report['missing_states'])} spectator states with coverage gaps")
    raise SystemExit(1 if report['missing_states'] else 0)
