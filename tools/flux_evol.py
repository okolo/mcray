import numpy as np
import sys
import os

import matplotlib.pyplot as plt

params = {'weight_col': 2, 'delay_col': 6, 'distance_Mpc': 1427.48, 'step_h': 0.5, 'max_steps': 100}

def usage_exit(reason=''):
    print(reason, file=sys.stderr)
    print("usage: python", sys.argv[0], 'data_file', *['{}={}'.format(k, v) for k, v in params.items()], file=sys.stderr)
    print("data must be sorted by column delay_col")
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
    if constr == bool:
        value = (value.lower() == "false" or value == "0")
        value = not value
    else:
        value = constr(value)
    params[kv[0]] = value


weight_col = params['weight_col']
col = params['delay_col']
step = params['step_h']
max_steps = params['max_steps']

data = np.loadtxt(path)

x_scale = params['distance_Mpc']*2.8590868063e+10  # convert to hours
col-=1
weight_col-=1

delay = data[:, col] * x_scale
delay = delay - np.min(delay)

if weight_col>=0:
    assert weight_col < data.shape[1]
    weights = data[:,weight_col]
else:
    weights = np.ones(delay.shape)

j=0
for i in range(max_steps):
    k = np.searchsorted(delay[j:], step*(i+1), side='right')
    print(step*(i+0.5), np.sum(weights[j:k]), k-j)
    j = k
