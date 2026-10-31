#!/usr/bin/env python3
import numpy as np
import sys,os,time
from tqdm import tqdm
import threading

def getDivB(fname):
	print("Reading...",file=sys.stderr)
	B = np.loadtxt(fname)
	print("Processing...",file=sys.stderr)
	N = round(len(B)**(1/3))
	Bx = B[:,3].reshape((N,N,N))
	By = B[:,4].reshape((N,N,N))
	Bz = B[:,5].reshape((N,N,N))

	#x_s = list(set(B[:,0]))
	dx = (max(B[:,0])-min(B[:,0]))/(N-1)
	print("dx =",dx,file=sys.stderr)

	Bx_dx = (Bx[1:]-Bx[:-1])/dx
	Bx_dx_y = (Bx[:,1:]+Bx[:,:-1])/2
	Bx_dx_y_z = (Bx[:,:,1:]+Bx[:,:,:-1])/2

	By_x = (By[1:]+By[:-1])/2
	By_x_dy = (By[:,1:]-By[:,:-1])/dx
	By_x_dy_z = (By[:,:,1:]+By[:,:,:-1])/2

	Bz_x = (Bz[1:]+Bz[:-1])/2
	Bz_x_y = (Bz[:,1:]+Bz[:,:-1])/2
	Bz_x_y_dz = (Bz[:,:,1:]-Bz[:,:,:-1])/dx

	divB = Bx_dx_y_z + By_x_dy_z + Bz_x_y_dz

	s = str(divB.min())+'\t'
	s += str(divB.max())+'\t'
	s += str(divB.mean())+'\t'
	s += str(divB.std())+'\t'
	s += str(abs(divB).min())+'\t'
	s += str(abs(divB).max())+'\t'
	s += str(abs(divB).mean())+'\t'
	s += str(abs(divB).std())
	return s

fl = open("divB_and_waves_fixed.dat","w")
s = "step\t"
s += "divB.min\t"
s += "divB.max\t"
s += "divB.mean\t"
s += "divB.std\t"
s += "abs(divB).min\t"
s += "abs(divB).max\t"
s += "abs(divB).mean\t"
s += "abs(divB).std"
print(s,file=fl)

thrs = []

def thrProcessing(step):
	fname = "/tmp/steps_"+str(step)+"_z0.002106_mfdump/magneticField.out"
	s = str(step)+"\t"+getDivB(fname)
	print(s,file=fl)
	fl.flush()
	if step==1.0001:
		sys.exit(0)

for step in tqdm([1.3, 1.2, 1.1, 1.08, 1.06, 1.04, 1.025, 1.01, 1.008, 1.006, 1.005, 1.001, 1.0008, 1.0005, 1.0002, 1.0001]):
	cmd = "./CRbeam -p 10 -z 0.002106 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth 2 -os -oo -o /tmp/steps_"+str(step)+"_z0.002106_mfdump -jet 3.14 -N 1 -rseed 42 --lEGMF 1 --outmagnetic --accuracy 1 -mfstep "+str(step)
	#os.system(cmd)
	#thrs.append( threading.Thread(target=thrProcessing, args = (step,)) )
	#thrs[-1].daemon = True
	#thrs[-1].start()
	thrProcessing(step)

time.sleep(300*60)
