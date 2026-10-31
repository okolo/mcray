#!/usr/bin/env python3
import os,sys,time,subprocess,gc
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

if 'MF' in os.environ:
	fname = os.environ['MF']
elif len(sys.argv)>1:
	fname = sys.argv[1]
else:
	print("Usage: MF=/tmp/MFsaveTest_FFT_seed37_FILLZERO/magneticField.out PYTHONSTARTUP=\""+sys.argv[0]+"\" python3")
	exit(1)

print("Load file:",fname)
table = np.loadtxt(fname)
print("Calculating...")

N = round(len(table)**(1/3))

Bx = table[:,3].reshape((N,N,N))
By = table[:,4].reshape((N,N,N))
Bz = table[:,5].reshape((N,N,N))

absB = (Bx**2+By**2+Bz**2)**0.5

print("mean |B|",absB.mean()*1e9,"nG")

plt.matshow(absB[ int(N/2) ])
plt.title("z = "+str(table[int(N/2),0])+" Mpc")
plt.show(block=False)

plt.matshow(absB[:,int(N/2),: ])
plt.title("y = "+str(table[int(N/2),1])+" Mpc")
plt.show(block=False)

plt.matshow(absB[:,:,int(N/2) ])
plt.title("z = "+str(table[int(N/2),2])+" Mpc")
plt.show(block=False)

def pltfree():
	'''
	Clear pictures for free RAM
	'''
	gc.collect()
	plt.cla()
	plt.clf()	# normally I use these lines to release the memory
	plt.close()


for i in range(absB.shape[2]):
	plt.matshow(absB[:,:,i])
	plt.title("z = "+str(table[i,2])+" Mpc")
	plt.savefig(f"/tmp/field{i:03d}.png")
	pltfree()