#!/net/minus/home/dolgikh/anaconda3/bin/python3

from FileHelpers import *
import psutil

def z2D(z):
	'''
	Convert redshift z to distance in kpc
	'''
	redshiftV=3e5*((1+z)**2-1)/((1+z)**2+1) # km/sec
	H = 70 # km/(sec*Mpc)
	return redshiftV/H * 1e3 # in Kpc

def D2z(D):
	'''
	Convert distance to recommended redshift z (more than need)
	'''
	z=0.0000587
	while(z2D(z)<D):z*=1.01
	if D<28000:
		return round(z*1.06,7)
	else:
		return round(z*1.15,7)

# Very good void
obsx=-2265.95573062/10000.
obsy= 8863.84070364/10000.
obsz= 4037.04998821/10000.

obsr=1000 # kpc
D = 35000 # kpc, start distance
step=1500 # kpc, stepsize
tagname="multiAnother2Void"

jet_theta,jet_phi=xyz2thetaPhi(obsx,obsy,obsz)
jet_angle=min(np.pi/2,2*np.arctan(obsr/D))

N=3*1000*1000

import os,time

def runCRbeam():
	global obsx,obsy,obsz,obsr,D,jet_theta,jet_phi,jet_angle,N
	z=D2z(obsr+D)
	nth=0
	for i in range(5):
		time.sleep(2)
		nth+=int((1-psutil.cpu_percent()/100)*psutil.cpu_count())-1
	nth=int(nth/5)
	if(nth<1):nth=1
	if(nth>60):print("Cluster is free, very nice!")
	cmd="./CRbeam -p 10 -z "+str(z)+" -E 1e19 -e 1e15 -mft -mf 1e-9 -nth "+str(nth)+" -o output/"+tagname+"_D"+str(round(D))+"_R"+str(obsr)+" -os"
	cmd+=" -jet "+str(jet_angle)+" -jett "+str(jet_theta)+" -jetp "+str(jet_phi)
	cmd+=" -obs 1 -obsr "+str(obsr)+" -obsx "+str(obsx*D)+" -obsy "+str(obsy*D)+" -obsz "+str(obsz*D)+" -N "+str(N)+" -rseed 42 --lEGMF 1"
	print("RUN COMMAND: "+cmd)
	os.system(cmd)
	return "output/"+tagname+"_D"+str(round(D))+"_R"+str(obsr)+"/z0/proton"

while(D<80500):
	outpfile=runCRbeam() # Run CRbeam
	outpfile=CRbeamFile(outpfile) # Read result
	D+=step # Update distance
	tmpx=outpfile['source_x'][outpfile.in_target]
	tmpy=outpfile['source_y'][outpfile.in_target]
	tmpz=outpfile['source_z'][outpfile.in_target]
	jet_angle=0.1+3*approxThetaParam(tmpx,tmpy,tmpz) # Update jet angle
	jet_theta,jet_phi=outpfile.target_direction # Update direction on target
	if(outpfile.in_target_num<200):
		print("Частиц слишком мало, увеличиваем число частиц")
		N*=2
#		if(N>6000000):
#			N=20000000
#			jet_angle=0.01+2.5*approxThetaParam(tmpx,tmpy,tmpz)
#			print("Упираюсь в предел. Останавливаюсь на N=6 000 000 и уменьшаю угол джета: jet_angle =",jet_angle)
	if(outpfile.in_target_num<3):
		print("Частиц неприлично мало. Выхожу...")
		exit(0)
#		jet_angle*=0.85
#		D-=step
#		os.system("rm -rf \"output/"+tagname+"_D"+str(round(D))+"_R"+str(obsr)+"/\"")
	else:
		print("Частиц попало достаточно продолжаем исполнение!")

print("Задача завершена.")
