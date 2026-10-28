#!/usr/bin/env python3
"""Validate completed p/n generator runs before using them as reference tables."""
import argparse
import hashlib
import json
import math
from pathlib import Path

SPECIES = {'e': 0, 'ae': 1, 'ph': 2, 'en': 3, 'mn': 4,
           'aen': 6, 'amn': 7, 'n': 9, 'p': 10}


def validate_run(directory, primary, density, trials, kmax):
    status = directory / 'exit-code.txt'
    if not status.is_file() or status.read_text().strip() != '0':
        raise ValueError(f'Generator run is not recorded as successful: {directory}')
    summaries = {}
    common_kmin, common_cells = None, None
    for species, secondary in SPECIES.items():
        path = directory / f'{primary}-{species}'
        with path.open() as stream:
            tokens = (token for line in stream for token in line.split())
            header = [next(tokens) for _ in range(5)]
            part, sec, kmin, bins, events = int(header[0]), int(header[1]), float(header[2]), int(header[3]), int(header[4])
            if (part, sec, bins, events) != (SPECIES[primary], secondary, density, trials):
                raise ValueError(f'Incompatible header: {path}')
            if not math.isfinite(kmin) or kmin <= 0 or kmin > kmax:
                raise ValueError(f'Invalid threshold: {path}')
            cells, total = 0, 0
            while True:
                try:
                    low = int(next(tokens))
                except StopIteration:
                    break
                high = int(next(tokens))
                cells += 1
                if low == -1 and high == -1:
                    continue
                if not 0 <= low <= high <= 8*density:
                    raise ValueError(f'Invalid secondary range: {path}, cell {cells}')
                for _ in range(high-low+1):
                    count = int(next(tokens))
                    if count < 0:
                        raise ValueError(f'Negative count: {path}')
                    total += count
        expected_cells = math.floor(math.log10(kmax/kmin)*density)+1
        if cells != expected_cells:
            raise ValueError(f'Incomplete photon grid: {path}: {cells} != {expected_cells}')
        if common_kmin is None:
            common_kmin, common_cells = kmin, cells
        elif (kmin, cells) != (common_kmin, common_cells):
            raise ValueError(f'Secondary files use different photon grids: {path}')
        summaries[path.name] = {'photon_cells': cells, 'total_generated_count': total,
                               'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    cross_section = directory / primary
    cells = 0
    with cross_section.open() as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(f'Invalid cross-section row: {cross_section}')
            energy, sigma = map(float, fields)
            expected_energy = common_kmin*10**(cells/density)
            if not math.isfinite(energy) or not math.isfinite(sigma) or sigma < 0 or abs(energy/expected_energy-1) > 1e-5:
                raise ValueError(f'Invalid cross-section grid: {cross_section}')
            cells += 1
    if cells != common_cells:
        raise ValueError(f'Incomplete cross-section file: {cross_section}')
    summaries[primary] = {'photon_cells': cells,
                         'sha256': hashlib.sha256(cross_section.read_bytes()).hexdigest()}
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p-directory', type=Path, required=True)
    parser.add_argument('--n-directory', type=Path, required=True)
    parser.add_argument('--bins-per-decade', type=int, required=True)
    parser.add_argument('--trials-per-cell', type=int, required=True)
    parser.add_argument('--kmax-gev', type=float, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    args = parser.parse_args()
    if args.bins_per_decade < 1 or args.trials_per_cell < 1 or not math.isfinite(args.kmax_gev) or args.kmax_gev <= 0:
        parser.error('Invalid generation parameters')
    try:
        result = {primary: validate_run(directory, primary, args.bins_per_decade,
                                       args.trials_per_cell, args.kmax_gev)
                  for primary, directory in [('p', args.p_directory), ('n', args.n_directory)]}
    except (OSError, ValueError, StopIteration) as error:
        parser.error(f'Table validation failed: {error}')
    # Publish a manifest only after both complete runs pass; do not label partial runs valid.
    with args.manifest.open('x') as output:
        json.dump({'bins_per_decade': args.bins_per_decade,
                   'trials_per_cell': args.trials_per_cell,
                   'kmax_GeV': args.kmax_gev, 'files': result}, output, indent=2)
        output.write('\n')
    print('Validated complete p/n runs: 18 secondary tables and 2 cross-section files')


if __name__ == '__main__':
    main()
