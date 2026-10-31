#!/usr/bin/env python3
import os,sys,time,subprocess
import numpy as np

# initial params
P=0.1 # as in paper
rseed=42 # random number
N=1000000 # all number of particles
batch_size=20000
a=np.arctan(0.01) # detector size
kappa=np.log(1.-P)/(a-1.) # initial kappa
theta=0 # initial angles, doesn't matter in first run
phi=0

def ttt():
	return time.strftime("[%H:%M:%S]: ")

def bashExec(q):
	return subprocess.check_output(q,shell=True).decode("UTF-8")

def xyz2thetaPhi(x,y,z):
	r=(x*x+y*y+z*z)**0.5
	theta=np.arccos(z/r)
	phi=1.5707963267948957
	try:
		phi=np.arctan(y/x)+(x<0)*3.1415926535897932
		phi+=(phi<0)*2.*3.1415926535897932
	except:
		print("Warning: singular phi")
	return theta,phi

isotropicProtons="./CRbeam -N 1000000 -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth 32 -obs 1 -obsr 100 -oo -o isotropicProtons -jet 3.14159265359"

def updateParams(fname):
	'''
	Update params by equations (3.5)-(3.6)
	'''
	prot=np.loadtxt(fname)
	x=prot[:,2] # Destination coords
	y=prot[:,3]
	z=1-prot[:,4]
	in_target=prot[:,9].astype(np.bool) # Is we in target?
	if(sum(in_target)==0):
		print("No particles in target. Return old values")
		global theta,phi,kappa
		return theta,phi,kappa
	weights=prot[:,1][in_target] # weights
	source_x=prot[:,10][in_target] # Start direction
	source_y=prot[:,11][in_target]
	source_z=prot[:,12][in_target]
	mu_x=sum(source_x*weights)
	mu_y=sum(source_y*weights)
	mu_z=sum(source_z*weights)
	mu_n=(mu_x**2+mu_y**2+mu_z**2)**0.5
	mu_x/=mu_n
	mu_y/=mu_n
	mu_z/=mu_n
	a=sum( weights*(mu_x*source_x+mu_y*source_y+mu_z*source_z)**2 )
	a/=sum(weights)
	newkappa=np.log(1.-P)/(a-1.)
	newtheta,newphi=xyz2thetaPhi(mu_x,mu_y,mu_z)
	print("Update params: (theta,phi)=("+str(newtheta)+","+str(newphi)+")")
	print("               a = "+str(a))
	print("               kappa = "+str(newkappa))
	return newtheta,newphi,newkappa

def weightf(kappa,scal):
	return kappa*np.exp(kappa*(scal-1))/(2*3.1415926535897932*(1-np.exp(-2*kappa)))

def runCRbeam(theta,phi,kappa,fid,cutoff_weight=0.01):
	maxweight=kappa/(2*3.1415926535897932*(1-np.exp(-2*kappa)))
	cutoff_weight=cutoff_weight*maxweight
	cutoffangle=3.1415
	for i in range(314):
		if(weightf(kappa,np.cos(float(i)/100.))<cutoff_weight):
			cutoffangle=float(i)/100.
			break
	if(cutoffangle<3.14):
		print("Cutoff angle:"+str(cutoffangle))
	s="./CRbeam -N "+str(batch_size)
	s+=" -p 10 -z 0.00234 -E 1e19 -e 1e15 -mft -mf 1e-9 -nth 32 -obs 1 -obsr 100 -os -o outputs/vmfProtons_"+fid
	s+=" -jet "+str(cutoffangle)
	s+=" -jett "+str(theta)
	s+=" -jetp "+str(phi)
	s+=" -jetprm "+str(kappa)
	s+=" -rseed "+str(rseed)
	print("Run CRbeam as:\n$ "+s)
	os.system(s)

iterations=int(N/batch_size)

for itr in range(iterations):
	fid=str(itr)
	ret=runCRbeam(theta,phi,kappa,fid,cutoff_weight=0.01)
	theta,phi,kappa=updateParams("outputs/vmfProtons_"+fid+"/z0/proton")
	os.system("cat outputs/vmfProtons_*/z0/proton >> outputs/allVMFprotons")
