import sys
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv)<2:
    print('usage: python', sys.argv[0],'w_dt_column_file [dt] [maxT]', file=sys.stderr)
    print('w_dt_column_file should contain weights in first column and delay in the second column', file=sys.stderr)
    print('w_dt_column_file should be ordered in dt: (sort -g -k 2)', file=sys.stderr)
    sys.exit(1)

data = np.loadtxt(sys.argv[1])
deltaT = data[:,1]
weight = data[:,0]

t_min = deltaT.min()
t_max = deltaT.max()
deltaT -= t_min
t_max -= t_min


dt=-10. # 10 bins by default

if len(sys.argv)>2:
    dt = float(sys.argv[2])
if dt<0:
    dt = (t_min-t_max)/dt

if len(sys.argv)>3:
    t_max = float(sys.argv[3])

n_bins = int((t_max-t_min)/dt + 0.5)
t_binning = np.linspace(t_min, t_max, n_bins+1)
f=np.zeros(n_bins)
j=0
for i, t in enumerate(t_binning[1:]):
    k=np.searchsorted(deltaT[j:], t, side='right')
    delT = (t-t_binning[i])
    f[i] = np.sum(weight[j:j + k])/delT
    j += k

t_binning = t_binning[:-1] + 0.5*dt # bin centers
plt.plot(t_binning, f)
plt.savefig(sys.argv[1] + ".pdf")
