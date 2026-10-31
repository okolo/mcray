import numpy as np
import sys
import os

import matplotlib.pyplot as plt

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

params = {'weight_col': 2, 'delay_col': 6, 'distance_Mpc': 1427.48, 'n_points': 40, 'logscale': False, 'suffix': '',
          'EmaxTeV': 1., 'psf': 0.4, 'right_percentile': 90., 'manual': False, 'cut_0' : False, 'bin_h' : 0.}

def usage_exit(reason=''):
    print(reason, file=sys.stderr)
    print("usage: python", sys.argv[0], 'data_file', *['{}={}'.format(k, v) for k, v in params.items()], file=sys.stderr)
    exit(1)

if len(sys.argv) < 2:
    usage_exit()

path = os.path.expanduser(sys.argv[1])


for par in sys.argv[2:]:
    kv=par.split('=')
    if len(kv)!=2:
        usage_exit('invalid argument ' + par)
    if kv[0] not in params:
        usage_exit('unknown parameter ' + kv[0])
    constr = type(params[kv[0]])
    value = kv[1]
    if constr == bool:
        value = (value.lower() == "false" or value == "0")
        value = not value
    else:
        value = constr(value)
    params[kv[0]] = value


weight_col = params['weight_col']-1
col = params['delay_col']-1
logscale = params['logscale']


data = np.loadtxt(path)

Emax = params['EmaxTeV']
data = data[data[:,0] <= Emax*1e12]
psf = params['psf']
data = data[data[:,2] <= psf]
cut0 = params['cut_0']
if cut0:
    data = data[data[:, col] != 0.]

rp = params['right_percentile'] # right percentile

suffix = params['suffix']

if params['manual']:
    suffix += "_manual"
if cut0:
    suffix += "cut0"

bin_h = params['bin_h']

if bin_h<=0.:
    out_name = '{0}_f{1}_Emax{2:g}_th{3}_p{4:g}{5}.png'.format(path, col+1, Emax, psf, rp , suffix)
else:
    out_name = '{0}_f{1}_Emax{2:g}_th{3}_b{4:g}{5}.png'.format(path, col + 1, Emax, psf, bin_h, suffix)

x_label = 't, [h]'
y_label = 'dN/dt [a.u.]'
x_scale = params['distance_Mpc']*2.8590868063e+10  # convert to hours


#### customize x_scale and label i.e.
#x_label = r'$\theta$'
#x_scale = 180.
#x_shift = -1.
# to use latex symbols uncomment the line below
#plt.rc('text', usetex=True)


delay = data[:, col] * x_scale
delay = delay - np.min(delay)
# the histogram of the data
plt.xlabel(x_label)
plt.ylabel(y_label)

# automatically select range
lp = 0  # left percentile

if logscale:
    xMin = np.min(delay[delay>0])
else:
    xMin = np.nanpercentile(delay, lp)

n_points = params['n_points']

if bin_h<=0.:
    xMax = np.nanpercentile(delay, rp)
else:
    xMax = xMin + n_points*bin_h

if logscale:
    bins = np.logspace(xMin, xMax, n_points)
else:
    bins = np.linspace(xMin, xMax, n_points)

if weight_col>=0:
    assert weight_col < data.shape[1]
    weights = data[:,weight_col]
else:
    weights = np.ones(bins.shape)

#average = np.sum(delay*weights)/np.sum(weights)
median = weighted_quantile(delay,[0.5],sample_weight=weights)[0]

if logscale:
    step = (xMax / xMin) ** (1. / (n_points - 1))
    lg_step = np.log(step)
    indexes = (np.log(delay) / lg_step).astype(np.int)
    w = np.power(1. / step, indexes)
    weights *= w

if not params['manual']:
    plt.hist(delay, bins, weights=weights, normed=1, facecolor='green', alpha=0.5, label='<dT>=' + '{0:g} h'.format(median))
else:
    dT1 = bins[0]
    gr_bins = np.zeros(2*len(bins))
    gr_bins[0] = 0.999*dT1
    gr_vals = np.zeros(2*len(bins))
    for i, dT2 in enumerate(bins[1:]):
        gr_bins[1 + 2 * i] = 1.001 * dT1
        gr_bins[2 + 2 * i] = 0.999 * dT2
        sub_data = data[dT1/x_scale <= data[:,col]]
        sub_data = sub_data[sub_data[:, col]<dT2/x_scale]
        w = np.sum(sub_data[:,weight_col])
        gr_vals[1 + 2 * i] = w
        gr_vals[2 + 2 * i] = w
        dT1 = dT2

    plt.plot(gr_bins, gr_vals, label='<dT>=' + '{0:g} h'.format(median))

plt.legend(loc='upper right')

plt.savefig(out_name)
plt.close()


