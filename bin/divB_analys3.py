import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
x = np.unique(data[:,0])
y = np.unique(data[:,1])
z = np.unique(data[:,2])
nx, ny, nz = len(x), len(y), len(z)

Bx = data[:,3].reshape(nx, ny, nz)
By = data[:,4].reshape(nx, ny, nz)
Bz = data[:,5].reshape(nx, ny, nz)

dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]

dBx_dx = np.gradient(Bx, dx, axis=0)
dBy_dy = np.gradient(By, dy, axis=1)
dBz_dz = np.gradient(Bz, dz, axis=2)

divB = dBx_dx + dBy_dy + dBz_dz
print("Статистика |div B|:", np.min(np.abs(divB)), np.mean(np.abs(divB)), np.max(np.abs(divB)))
