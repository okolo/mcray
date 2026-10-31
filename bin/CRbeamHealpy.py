#!/usr/bin/env python3
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys

ENABLE_FLIP=False#True
rot=(0,90,0) #
enableLog=True #
enableRef=True #
enableGrad=True #
enableBigTargetFix=True #
nside=512

fname="photon"
try:
	fname=sys.argv[1]
except:
	print("No file specified, try to open \"photon\"")

print("Loading file...")
Norm=1.
with open(fname,"r") as fl: #
	tmp=fl.readline()
	import re
	Norm=float(re.compile("T/Mpc = ([0-9\.]+)").findall(tmp)[0])*1e3

phot=np.loadtxt(fname)
x=phot[:,2]*Norm
y=phot[:,3]*Norm
z=(1-phot[:,4])*Norm
print("Max distance (by datafile)",max((x*x+y*y+z*z)**0.5))
Px=phot[:,5]
Py=phot[:,6]
Pz=1-phot[:,7]
weights=phot[:,1]
print("Min weight:",weights.min(),"Max weight:",weights.max())
if(weights.max()-weights.min()>1e+2):
	weights=np.log(1+weights)
	print("Log weights for plot")
elif(weights.max()-weights.min()>1e-6):weights-=weights.min()
weights/=weights.max()
#weights=np.ones_like(weights)
print("Update: Min weight:",weights.min(),"Max weight:",weights.max())
in_target=np.zeros_like(weights).astype(np.bool)
try:
	in_target=phot[:,9].astype(np.bool)
	print("All particles:",len(in_target),", in target:",sum(in_target))
except:
	print("All particles:",len(in_target))
	print("No target column found")


params_dict=dict()
maxDistance=0
try:
	import pathlib
	p=pathlib.Path(fname).resolve().parent.parent.as_posix()
	p+="/params"
	params=open(p.replace("//","/"),'r').readlines()
	print("Load run params from",p,"file")
	for l in params:
		if(l.startswith("#")):continue
		if(l.startswith("--")):
			par=l.replace("\n","").split("\t")
			try:
				params_dict[par[0][2:]]=par[1]
			except IndexError:
				params_dict[par[0][2:]]=""
	params_dict['z']=float(params_dict['z'])
	if(params_dict['z']>1):
		print("Warning: redshift is so big! You cannot use z->D approximation!")
	redshiftV=3e5*((1+params_dict['z'])**2-1)/((1+params_dict['z'])**2+1) # km/sec
	H = 70 # km/(sec*Mpc)
	maxDistance=redshiftV/H * 1e3 # in Kpc
	print("Max distance (by redshift):",maxDistance,"Kpc")
except Exception as e:
	print("Warning: failed to parse run commands:",e)


target_x=None
target_y=None
target_z=None
#target_r=None
if("observer-x" in params_dict and sum(in_target)>0):
	target_x=float(params_dict["observer-x"])
	target_y=float(params_dict["observer-y"])
	target_z=float(params_dict["observer-z"])
#	target_r=float(params_dict["observer-radius"])
	print("Automatically found target coords:",(target_x,target_y,target_z))
#	target_x/=maxDistance#/Norm
#	target_y/=maxDistance#/Norm
#	target_z/=maxDistance#/Norm
#	target_r/=maxDistance#/Norm
#	print("Normilized on",maxDistance,"coords:",(target_x,target_y,target_z))


lbl=""
try:
	lbl=sys.argv[1][:sys.argv[1].index("/")]
except:
	lbl=sys.argv[1]

def xyz2thetaPhi(x,y,z):
	if(ENABLE_FLIP):
		x+=z
		z=x-z
		x-=z
		z*=-1
	r=(x*x+y*y+z*z)**0.5
	theta=np.arccos(z/r)
	phi=1.5707963267948957
	try:
		phi=np.arctan(y/x)+(x<0)*3.1415926535897932
		phi+=(phi<0)*2.*3.1415926535897932
	except:
		print("Warning: singular phi")
	return theta,phi

def calcTargetCoords(x_in,y_in,z_in,verbose=True):
	global target_x,target_y,target_z
	r=max(x_in*x_in+y_in*y_in+z_in*z_in)**0.5
	x0=np.array([x_in.mean(),y_in.mean(),z_in.mean()])
	x0=x0/np.linalg.norm(x)*r
	if(not target_x is None):
		trueRes=np.array([target_x,target_y,target_z])
		if( sum((x0-trueRes)**2)**0.5>r ):
			print("Warning: target coordinats inconsistent. Check params file!")
			print("Warning: by file:",[target_x,target_y,target_z])
			print("Warning: by calc:",list(x0))
		return trueRes
	return x0
	x0dev=max(x_in.std(),y_in.std(),z_in.std())
	foo=lambda x: ((x_in-x[0])**2+(y_in-x[1])**2+(z_in-x[2])**2).std()/r+( (x-x0).dot(x-x0)/x0dev**2 )
	if(verbose):print("Try to found target coords. Start coords:",list(x0),"Err:",foo(x0))
	from scipy.optimize import minimize
	res=minimize(foo,x0);
	if(verbose):print(res)
	return res.x

source_x=None
source_y=None
source_z=None

try:
	source_x=phot[:,10]*Norm
	source_y=phot[:,11]*Norm
	source_z=phot[:,12]*Norm
	theta,phi=xyz2thetaPhi(source_x.mean(),source_y.mean(),source_z.mean())
	print("Jet direction: theta =",theta,", phi =",phi)
	theta,phi=xyz2thetaPhi(source_x[in_target].mean(),source_y[in_target].mean(),source_z[in_target].mean())
	print("Direction on target: theta =",theta,", phi =",phi)
except:
	print("Source direction are not specified")


import healpy as hp
from tqdm import tqdm
print("nside =",nside)
npix = hp.nside2npix(nside)
# Go from HEALPix coordinates to indices
thetas,phis=xyz2thetaPhi(x,y,z)
indices = hp.ang2pix(nside, thetas, phis)

# Initate the map and fill it with the values
hpxmap = np.zeros(npix, dtype=np.float)
tmp=weights
for i in tqdm(range(len(tmp))):
    hpxmap[indices[i]] += tmp[i]

if(enableLog):hpxmap=np.log(1.+hpxmap)

# Inspect the map
hp.mollview(hpxmap,rot=rot)
if(enableGrad):hp.graticule()
if(enableRef):
	t,p=xyz2thetaPhi(0., 0., 1.)
	if(sum(in_target)>0):
		target_coords=calcTargetCoords(x[in_target],y[in_target],z[in_target],verbose=False)
		t,p=xyz2thetaPhi(target_coords[0],target_coords[1],target_coords[2])
		if(np.isnan(p) or np.isinf(p)):p=0.0
	print("Mark target(if exist) as red point (theta,phi) =",(t,p))
	hp.projscatter(t, p,s=70,color='red',marker='.')
plt.title("End particles sphere by source "+lbl)
#plt.savefig(fname="/tmp/MyFigure"+str(nside)+".png",dpi=240,orientation='landscape')
plt.show(block=not (sum(in_target)>0))

if(sum(in_target)>0):
	#target_coords=[0,0,max(x[in_target]*x[in_target]+y[in_target]*y[in_target]+z[in_target]*z[in_target])**0.5] # TODO
	target_coords=calcTargetCoords(x[in_target],y[in_target],z[in_target])
	print("target_coords =",list(target_coords))
	npix = hp.nside2npix(nside)
	# Go from HEALPix coordinates to indices
	thetas,phis=xyz2thetaPhi(target_coords[0]-x[in_target],target_coords[1]-y[in_target],target_coords[2]-z[in_target])
	indices = hp.ang2pix(nside, thetas, phis)
	hpxmap = np.zeros(npix, dtype=np.float)
	tmp=weights[in_target]
	for i in tqdm(range(len(tmp))):
		hpxmap[indices[i]] += tmp[i]
	if(enableLog):hpxmap=np.log(1.+hpxmap)
	hp.mollview(hpxmap,rot=rot)
	if(enableGrad):hp.graticule()
	if(enableRef):hp.projscatter(t, p,s=70,color='red',marker='.')
	plt.title("Where on the target observer see "+lbl)
	plt.show(block=False) ######### PLOT
	# Source directions
	npix = hp.nside2npix(nside)
	# Go from HEALPix coordinates to indices
	thetas,phis=xyz2thetaPhi(source_x,source_y,source_z)
	indices = hp.ang2pix(nside, thetas, phis)
	hpxmap = np.zeros(npix, dtype=np.float)
	tmp=in_target+1
	for i in tqdm(range(len(tmp))):
		hpxmap[indices[i]] = max(tmp[i],hpxmap[indices[i]])
	hp.mollview(hpxmap,rot=rot)
	if(enableGrad):hp.graticule()
	if(enableRef):hp.projscatter(t, p,s=70,color='red',marker='.')
	plt.title("Initial particles directions "+lbl)
	plt.show(block=False) ######### PLOT
	# Direction
	npix = hp.nside2npix(nside)
	# Go from HEALPix coordinates to indices
	thetas,phis=xyz2thetaPhi(Px[in_target],Py[in_target],Pz[in_target])
	if(enableBigTargetFix):
		print("Fix directions...")
		thetasFix,phisFix=xyz2thetaPhi(x[in_target],y[in_target],z[in_target])
		thetas-=thetasFix
		phis-=3.1415926535897932*(thetas<0)
		phis+=2*3.1415926535897932*(phis<0)
		thetas=np.abs(thetas)
		#phis-=phisFix
		print("Fixed thetas: ",thetas.min(),thetas.max())
	indices = hp.ang2pix(nside, thetas, phis)
	print(indices)
	hpxmap = np.zeros(npix, dtype=np.float)
	tmp=weights[in_target]
	for i in tqdm(range(len(tmp))):
		hpxmap[indices[i]] += tmp[i]
	if(enableLog):hpxmap=np.log(1.+hpxmap)
	hp.mollview(hpxmap,rot=rot)
	if(enableGrad):hp.graticule()
	if(enableRef):hp.projscatter(0.0, xyz2thetaPhi(0.,0.,1.)[1],s=70,color='red',marker='.')
	plt.title("What observer see "+lbl)
	plt.show(block=True)#'''
