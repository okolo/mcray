#!/usr/bin/env python3
import os,sys,time,subprocess,gc
import numpy as np
from tqdm import tqdm

fname = "output/MFsaveTest_FFT_seed37_v2/magneticField.out"
#fname = "output/MFsaveTest/magneticField.out"
#fname = "output/MFsaveTest_FFT/magneticField.out"

if 'MF' in os.environ:
	fname = os.environ['MF']

if len(sys.argv)>1:
	fname = sys.argv[1]

print("Load file:",fname)
table = np.loadtxt(fname)
print("Calculating...")

N = round(len(table)**(1/3))

Bx = table[:,3].reshape((N,N,N))
By = table[:,4].reshape((N,N,N))
Bz = table[:,5].reshape((N,N,N))

absB = (Bx**2+By**2+Bz**2)**0.5

print("mean |B|",absB.mean()*1e9,"nG")
print("Calculating FFT...")

Bx_k = np.fft.rfftn(Bx,axes=(0,1,2))
By_k = np.fft.rfftn(By,axes=(0,1,2))
Bz_k = np.fft.rfftn(Bz,axes=(0,1,2))

Bx_k_abs = (abs(Bx_k)**2+abs(By_k)**2+abs(Bz_k)**2)**0.5

print("approx mean |B| by image",((Bx_k_abs**2).mean()/N**3)**0.5 *1e9,"nG")

# By DeepSeek

weights_z = np.ones(N//2+1, dtype=float)
if N % 2 == 0:
    weights_z[0] = 1.0
    weights_z[-1] = 1.0
    weights_z[1:-1] = 2.0
else:
    weights_z[0] = 1.0
    weights_z[1:] = 2.0

# Преобразуем в трёхмерный массив весов
weights = weights_z.reshape(1, 1, -1)  # (1,1,N//2+1)

# Суммируем с весами
E_fourier = np.sum(weights * (np.abs(Bx_k)**2 + np.abs(By_k)**2 + np.abs(Bz_k)**2)) / N**3

print("true mean |B| by image",(E_fourier/N**3)**0.5 *1e9,"nG")


dx = table[1,2]-table[0,2] # Mpc
# Частоты
kx = 2*np.pi * np.fft.fftfreq(N, d=dx)
ky = 2*np.pi * np.fft.fftfreq(N, d=dx)
kz = 2*np.pi * np.fft.rfftfreq(N, d=dx)

# Сетка k-векторов
Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')
K = np.sqrt(Kx**2 + Ky**2 + Kz**2)

#i,j,kk = 20,20,20
#B_v = np.array([Bx_k[i,j,kk],By_k[i,j,kk],Bz_k[i,j,kk]])
#k_v = np.array([kx[i],ky[j],kz[kk]])

# Энергия в k-пространстве (без весов)
Pk = np.abs(Bx_k)**2 + np.abs(By_k)**2 + np.abs(Bz_k)**2

k_flat = K.flatten()
Pk_flat = Pk.flatten()

# Определяем бины, например, по логарифмической шкале от k_min до k_max
k_min = 2*np.pi / (N*dx)   # минимальная ненулевая частота
k_max = np.sqrt(3)*(np.pi/dx)  # максимальная частота (по диагонали)
bins = np.logspace(np.log10(k_min), np.log10(k_max), 200)
k_centers = (bins[1:] + bins[:-1]) / 2

print("Preparing spectrum")
Pk_binned = []
Pk_std_binned = []
for i in tqdm(range(len(bins)-1)):
	mask = (k_flat >= bins[i]) & (k_flat < bins[i+1])
	if np.sum(mask) > 0:
		Pk_binned.append(np.mean(Pk_flat[mask]))
		Pk_std_binned.append(np.std(Pk_flat[mask]))
	else:
		Pk_binned.append(0.0)
		Pk_std_binned.append(0.0)

Pk_binned = np.array(Pk_binned)
Pk_std_binned = np.array(Pk_std_binned)

def getDiv(Bx,By,Bz):
	dBx     = (   Bx[:-1,:,:] -  Bx[1:,:,:])/dx
	dBx_y   = (  dBx[:,:-1,:] + dBx[:,1:,:])/2
	dBx_y_z = (dBx_y[:,:,:-1]+dBx_y[:,:,1:])/2
	dBy     = (   By[:,:-1,:] -  By[:,1:,:])/dx
	dBy_x   = (  dBy[:-1,:,:] + dBy[1:,:,:])/2
	dBy_x_z = (dBy_x[:,:,:-1]+dBy_x[:,:,1:])/2
	dBz     = (   Bz[:,:,:-1] -  Bz[:,:,1:])/dx
	dBz_x   = (  dBz[:-1,:,:] + dBz[1:,:,:])/2
	dBz_x_y = (dBz_x[:,:-1,:]+dBz_x[:,1:,:])/2
	divB = dBx_y_z+dBy_x_z+dBz_x_y
	return divB

def getDiv2(Bx,By,Bz):
	dBx_dx = np.gradient(Bx, dx, axis=0)
	dBy_dy = np.gradient(By, dx, axis=1)
	dBz_dz = np.gradient(Bz, dx, axis=2)
	divB = dBx_dx + dBy_dy + dBz_dz
	return divB

print("Calculating divergence...")
#divB = getDiv(Bx,By,Bz)
#divB = getDiv(Bz,By,Bx)
#divB = getDiv(By,Bz,Bx)
divB = getDiv(Bz,Bx,By)
abs_divB = abs(divB)
print("min|divB| =",abs_divB.min())
print("mean|divB| =",abs_divB.mean())
print("max|divB| =",abs_divB.max())

import matplotlib.pyplot as plt
import time

plt.loglog(k_centers, Pk_binned, 'o-')
plt.loglog(k_centers, Pk_binned-Pk_std_binned, '.')
plt.loglog(k_centers, Pk_binned+Pk_std_binned, '.')
plt.loglog(k_centers, Pk_binned.max()*k_centers**(-11/3), '-')

plt.xlabel('k [1/Mpc]')
plt.ylabel('P(k)')  # размерность зависит от нормировки
plt.grid(True)
plt.show()

print("K.shape =",K.shape)
print("kz =",kz)

plt.matshow(np.log(Bx_k_abs[:,:,0]))
plt.title("kz = 0")
plt.show(block = False)

plt.matshow(np.log(Bx_k_abs[:,:,10]))
plt.title("kz = 10")
plt.show(block = False)


def pltfree():
	'''
	Clear pictures for free RAM
	'''
	gc.collect()
	plt.cla()
	plt.clf()	# normally I use these lines to release the memory
	plt.close()


for i in range(K.shape[2]):
#	plt.matshow(K[:,:,i])
	plt.matshow(np.log(Bx_k_abs[:,:,i]))
	plt.title("kz = "+str(i))
#	plt.show(block=False)
	plt.savefig(f"/tmp/spectrum{i:03d}.png")
	pltfree()

for i in range(absB.shape[2]):
#	plt.matshow(K[:,:,i])
	plt.matshow(absB[:,:,i])
	plt.title("z = "+str(table[i,2])+" Mpc")
#	plt.show(block=False)
	plt.savefig(f"/tmp/field{i:03d}.png")
	pltfree()