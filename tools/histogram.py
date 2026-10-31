import numpy as np
import sys
import os

import matplotlib.pyplot as plt

weight_col=2
if len(sys.argv) < 3:
    print("usage: python", sys.argv[0], "file col [weight_col]\ndefault weight_col=" + str(weight_col), file=sys.stderr)
    exit(0)

path = os.path.expanduser(sys.argv[1])
col = int(sys.argv[2])
if len(sys.argv) > 3:
    weight_col = int(sys.argv[3])

data = np.loadtxt(path)
out_name = path + "_" + str(col) + ".png"

x_label = str(col)
x_scale = 4.0812892343e+13 # 1427.48 Mpc/hour
x_shift = 0.
col-=1
weight_col-=1

#### customize x_scale and label i.e.
#x_label = r'$\theta$'
#x_scale = 180.
#x_shift = -1.
# to use latex symbols uncomment the line below
#plt.rc('text', usetex=True)


real = data[:,col]*x_scale + x_shift
# the histogram of the data
plt.xlabel(x_label)

# automatically select range
lp = 5 # percentile
rp = 40
xMin = np.nanpercentile(real,lp)
xMax = np.nanpercentile(real,100-rp)

bins = np.linspace(xMin, xMax, 100)
weights=None
if weight_col>=0:
    assert weight_col < data.shape[1]
    weights = data[:,weight_col]

plt.hist(real, bins, weights=weights, normed=1, facecolor='green', alpha=0.5, label=str(col+1))

plt.legend(loc='upper right')

plt.savefig(out_name)
plt.close()


