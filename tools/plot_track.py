import sys
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

legend = {0: 'electron', 1: 'positron', 2: 'photon', 3: 'proton'}

def plot_colourline(x,y,fc):
    fc[fc==10]=3 # make more contrast color scale
    c = cm.jet((fc-np.min(fc))/(np.max(fc)-np.min(fc)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        label = legend.pop(fc[i], None)
        if label is None:
            ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])
        else:
            ax.plot([x[i], x[i + 1]], [y[i], y[i + 1]], c=c[i], label=label)
    return

if len(sys.argv)<2:
    print('usage: python', sys.argv[0],'track_log', file=sys.stderr)
    print('track_log is output produced by parse_track.py', file=sys.stderr)
    sys.exit(1)

data = np.loadtxt(sys.argv[1])

E=0.1*np.log10(data[0,5]/data[:,5])
t=data[:,3]
x=data[:,8] * t
y=data[:,9] * t
cosTh=(1.-data[:,10]) # z/t
cosTh[0] = 1 # fix artifact
z=cosTh*t
dt=t-np.sqrt(x*x+y*y+z*z)
part_type=data[:,2].astype(np.int)

Th=np.arccos(cosTh)*180/np.pi

plot_colourline(t, Th, part_type)
plot_colourline(t, E, part_type)
#plt.plot(t, np.sqrt(x*x+y*y), 'xy')
# plt.plot(t, x, c=part_type, label='x')
# plt.plot(t, y, c=part_type, label='y')
#plt.plot(t, z, label='z')
#plt.plot(t, dt, label='dt')
plt.legend(loc='upper right')
plt.savefig(sys.argv[1] + ".pdf")
