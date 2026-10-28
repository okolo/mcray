"""Genealogical pilot statistics. Never issues physical acceptance."""
from collections import Counter
import math
import statistics


def summarize(path, job, protocol):
    histories, metadata = {}, {}
    complete = False
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            metadata.update(t.split('=', 1) for t in line.split() if '=' in t)
            if line.startswith('# batch=0 primaries='):
                complete = int(line.split('primaries=')[1]) == job['histories']
            continue
        if not line:
            continue
        batch, h, species, weight, energy = line.split()
        h, weight, energy = int(h), float(weight), float(energy)
        if (batch != '0' or not 0 <= h < job['histories'] or
                not math.isfinite(weight) or weight <= 0 or
                not math.isfinite(energy) or energy < 0):
            raise ValueError('Invalid weighted product row')
        row = histories.setdefault(h, Counter())
        row['yield:'+species] += weight
        row['energy:'+species] += energy  # Stored energy already includes weight.
        if species.startswith('A'):
            a, z = map(int, species[1:].split('Z'))
            if a < 2 or not 0 <= z <= a:
                raise ValueError('Invalid isotope')
            row['bound_count'] += weight
            row['bound_A_sum'] += a*weight
            row['bound_Z_sum'] += z*weight
            row['energy:bound'] += energy
            row['survival'] += weight*(a == 56 and z == 26)
        else:
            a, z = {'proton': (1, 1), 'neutron': (1, 0), 'electron': (0, -1),
                    'positron': (0, 1), 'neutrinoE': (0, 0), 'neutrinoAE': (0, 0),
                    'photon': (0, 0)}[species]
        row['baryons'] += a*weight
        row['charge'] += z*weight
        if species == 'photon':
            unweighted_energy = energy/weight
            if not math.isfinite(unweighted_energy):
                raise ValueError('Nonfinite unweighted photon energy')
            k = next(i for i in range(4) if protocol['photon_bin_edges_eV'][i+1] is None or
                     unweighted_energy < protocol['photon_bin_edges_eV'][i+1])
            row[f'photon_bin_{k}_count'] += weight
            row[f'photon_bin_{k}_energy_eV'] += energy
    if (not complete or set(histories) != set(range(job['histories'])) or
            metadata.get('seed') != str(job['seed']) or
            float(metadata.get('L_mpc', 'nan')) != protocol['distance_mpc'] or
            metadata.get('PD_photons') not in ['1', 'On'] or
            metadata.get('decay_photons') not in ['1', 'On']):
        raise ValueError('Incomplete histories or changed pilot scope')
    primary = metadata.get('primary', 'A56Z26').rstrip(';')
    if primary != 'A56Z26' or metadata.get('process') not in [None, 'talys-decay']:
        raise ValueError('Changed initial isotope or process')
    if job['profile'] != 'ordinary':
        if (float(metadata.get('fraction', 'nan')) != job['fraction'] or
                int(metadata.get('copies', '-1')) != job['copies']):
            raise ValueError('Changed weighted proposal')
        if int(metadata.get('extra_continuations', '-1')) < 0:
            raise ValueError('Missing clone accounting')
        if int(metadata['extra_continuations']) > job['histories']*(job['copies']-1):
            raise ValueError('More than one root PD split per lineage')
    if job['code'] == 'mcray':
        if 'E_per_nucleon_eV' in metadata and float(metadata['E_per_nucleon_eV']) != protocol['energy_per_nucleon_eV']:
            raise ValueError('Changed initial energy')
        control = protocol['native_controls'][job['level']]
        step = metadata.get('background_step', metadata.get('background_step_decades', 'nan'))
        eps = metadata.get('eps', metadata.get('accuracy', 'nan'))
        if float(step) != control['step'] or float(eps) != control['eps']:
            raise ValueError('Changed native numerical settings')
    elif float(metadata.get('step_fraction', 'nan')) != protocol['reference_step_fractions'][job['level']]:
        raise ValueError('Changed reference numerical settings')
    keys = [*protocol['targets']['cmb/mcray']['means'], 'baryons', 'charge']
    means, variances, maxima, nonzero = {}, {}, {}, {}
    for key in keys:
        values = [histories[h][key] for h in range(job['histories'])]
        means[key] = statistics.mean(values)
        variances[key] = statistics.variance(values)
        maxima[key] = max(values)
        nonzero[key] = sum(v != 0 for v in values)
    return dict(means=means, history_variances=variances,
                max_history_contribution=maxima, nonzero_histories=nonzero,
                extra_continuations=int(metadata.get('extra_continuations', '0')))


def batch_statistics(batches):
    if len(batches) < 2:
        return dict(batches=len(batches), status='insufficient pilot batches')
    keys = batches[0]['summary']['means']
    result = {}
    for key in keys:
        means = [b['summary']['means'][key] for b in batches]
        se = statistics.stdev(means)/math.sqrt(len(means))
        result[key] = dict(mean=statistics.mean(means), standard_error=se,
            batch_variance=statistics.variance(means),
            zero_variance_unresolved=se == 0,
            nonzero_histories=sum(b['summary']['nonzero_histories'][key] for b in batches),
            maximum_history_contribution=max(b['summary']['max_history_contribution'][key] for b in batches))
    seconds = None if any(b.get('recovered_complete_output') for b in batches) else statistics.mean(b['seconds'] for b in batches)
    return dict(batches=len(batches), mean_batch_seconds=seconds,
                observables=result, status='pilot diagnostics only; no acceptance')
