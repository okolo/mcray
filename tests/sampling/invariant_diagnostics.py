"""Read deterministic s-quantile controls and repeated full-cost timings."""
import csv
import math
import statistics
from pathlib import Path

PROBABILITIES = (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, .001, .01, .1, .25, .5,
                 .75, .9, .99, .999, 1-1e-4, 1-1e-6, 1-1e-8, 1-1e-10, 1-1e-12)


def summarize(paths):
    cases = {}
    for source in paths:
        source = Path(source)
        lines = source.read_text().splitlines()
        header = dict(item.split('=', 1) for item in lines[0].removeprefix('# ').split())
        rows = list(csv.DictReader((line for line in lines if not line.startswith('#')), delimiter='\t'))
        if tuple(float(row['u']) for row in rows) != PROBABILITIES:
            raise ValueError('Incomplete deterministic quantile family: ' + str(source))
        numbers = [float(row[key]) for row in rows for key in
                   ('x', 'small_tail_reference', 'relative_probability_error', 'reference_change_relative')]
        numbers += [float(value) for value in header.values()]
        if not all(math.isfinite(value) for value in numbers):
            raise ValueError('Nonfinite deterministic diagnostic')
        if any(float(row['x']) <= 0 or float(row['small_tail_reference']) <= 0 for row in rows):
            raise ValueError('Invalid quantile or tail probability')
        if any(float(a['x']) >= float(b['x']) for a, b in zip(rows, rows[1:])):
            raise ValueError('Non-increasing quantiles')
        errors = [abs(float(row['relative_probability_error'])) for row in rows]
        changes = [abs(float(row['reference_change_relative'])) for row in rows]
        norm_error = abs(float(header['total']) / float(header['reference']) - 1)
        enclosure = next((dict(item.split('=') for item in line[len('# envelope '):].split())
                          for line in lines if line.startswith('# envelope ')), None)
        envelope_ok = True
        if enclosure is not None:
            enclosure = {key: float(value) for key,value in enclosure.items()}
            if not all(map(math.isfinite,enclosure.values())):
                raise ValueError('Nonfinite envelope diagnostic')
            envelope_ok = (0 < enclosure['lower'] <= enclosure['upper'] <= 1.010000001*enclosure['lower']
                           and 0 <= enclosure['max_ratio'] <= 1 and 0 <= enclosure.get('max_local_ratio',0) <= 1
                           and enclosure['checks'] > 0 and '# bounded_tail_proposals 19' in lines)
        timings = {mode: [] for mode in (('cdf', 'rejection', 'nr', 'bounded') if enclosure is not None else ('cdf', 'rejection', 'nr'))}
        repeats = {mode: set() for mode in timings}
        for line in lines:
            if line.startswith('# benchmark '):
                _, _, mode, repeat, count, seconds, checksum = line.split()
                if mode not in timings:
                    raise ValueError('Changed timing family')
                if int(count) != 128 or int(repeat) not in (0, 1, 2):
                    raise ValueError('Changed timing protocol')
                timings[mode].append(float(seconds) / int(count))
                repeats[mode].add(int(repeat))
        if any(len(values) != 3 or min(values) <= 0 or not all(map(math.isfinite, values))
               or repeats[mode] != {0, 1, 2} for mode, values in timings.items()):
            raise ValueError('Incomplete timing blocks')
        medians = {mode: statistics.median(values) for mode, values in timings.items()}
        # Numerical diagnostics, separate from the unchanged 1%/10%, alpha=1%
        # Monte Carlo acceptance family. Coarse/fine differences are empirical.
        reference_ok = max(changes) <= 1e-4
        accuracy_ok = max(errors) <= 1e-3 and norm_error <= 1e-4 and envelope_ok
        cases[source.stem] = {
            'status': 'within_control' if reference_ok and accuracy_ok else 'inconclusive_reference' if not reference_ok else 'deviation',
            'quantiles': len(rows), 'max_relative_small_tail_probability_error': max(errors),
            'max_relative_reference_change': max(changes), 'relative_normalization_error': norm_error,
            'seconds_per_draw_blocks': timings, 'median_seconds_per_draw': medians,
            'cdf_over_rejection_time': medians['cdf'] / medians['rejection'],
            'cdf_over_nr_time': medians['cdf'] / medians['nr'],
        }
        if enclosure is not None:
            cases[source.stem].update(envelope=enclosure,
                bounded_over_rejection_time=medians['bounded']/medians['rejection'],
                cdf_over_bounded_time=medians['cdf']/medians['bounded'])
    if not cases:
        raise ValueError('No invariant diagnostics')
    return {
        'scope': 'Numerical quantile controls and full per-draw construction cost; not a Monte Carlo acceptance result',
        'small_tail_relative_control': 1e-3, 'reference_change_guard': 1e-4,
        'status': 'within_control' if all(c['status'] == 'within_control' for c in cases.values()) else 'unresolved',
        'cases': cases,
    }
