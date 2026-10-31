#!/usr/bin/python3
#!/usr/bin/env python3
# Requires:
# sudo apt install python3-pip python3-numpy python3-matplotlib
# Anaconda is not recommended
import re,time

try:
	import pathlib
except ImportError:
	import pathlib2 as pathlib

import numpy as np
from tqdm import tqdm
import time
import healpy as hp

def fixBitFile(fname):
	'''
	Rewrite bit files
	'''
	outfl=fname+"_fixed"
	outfl=open(outfl,"w")
	with open(fname,"r") as fl:
		columns=0
		for il in fl:
			if("#" in il):c=outfl.write(il)
			if(columns==0):
				c=outfl.write(il)
				columns=il.count("\t")+il.count(" ")
			else:
				tmp=il.count("\t")+il.count(" ")
				if(tmp!=columns):
					#print("Bad line:",il)
					continue
				else:c=outfl.write(il)
	outfl.close()
	return fname+"_fixed"

def xyz2thetaPhi(x,y,z):
	'''
	Convert Cartesian coords to angles on sphere in radians
	Input: numbers or 1D NumPy arrays with the same length
	'''
	r=(x*x+y*y+z*z)**0.5
	theta=np.arccos(z/r)
	phi=1.5707963267948957
	try:
		phi=np.arctan(y/x)+(x<0)*3.1415926535897932
		phi+=(phi<0)*2.*3.1415926535897932
	except:
		print("Warning: singular phi")
	return theta,phi

def approxThetaParam(x,y,z):
	'''
	Find approximation theta param for jet's width for aiming to the target
	Input: numbers or 1D NumPy arrays with the same length
	'''
	r=(x*x+y*y+z*z)**0.5
	R=r.min()
	dx=x.max()-x.min()
	dy=y.max()-y.min()
	dz=z.max()-z.min()
	dr=max(dx,dy,dz)
	return np.arctan(dr/R)

def thetaPhi2hpmap(thetas, phis, nside = 256, weights = None):
	npix = hp.nside2npix(nside)
	# Go from HEALPix coordinates to indices
	indices = hp.ang2pix(nside, thetas, phis)
	# Initate the map and fill it with the values
	hpxmap = np.zeros(npix, dtype=np.float)
	if weights is None:
		weights = np.ones_like(thetas)
	for i in tqdm(range(len(weights))):
		hpxmap[indices[i]] += weights[i]
	return hpxmap

def plt_hpmap(hpxmap):
	xsize=max(round(1.05*(len(hpxmap)*2/3.14159)**0.5)*2,800)
	hp.mollview(hpxmap,rot=(0,90,0),norm='log',xsize=xsize)
	hp.graticule()

def approxThetaParam2(x,y,z):
	'''
	Find approximation theta param for jet's width for aiming to the target
	Input: numbers or 1D NumPy arrays with the same length
	More careful update for approxThetaParam
	'''
	vec=np.vstack([x,y,z]).T
	mean_vec=vec.mean(axis=0)
	mean_vec/=np.linalg.norm(mean_vec)
	cos_s=vec.dot(mean_vec)/np.linalg.norm(vec,axis=1)
	mincos=cos_s.min()
	return np.arccos(mincos)


class CRbeamFile:
	def __init__(self,fname,verbose=True,readAnyway=True,READTIMEOUT=500,project=False):
		'''
		Dataloader for CRbeam files
		Call example: f=CRbeamFile("z0/proton")
		verbose -- show debug messages in stdout
		readAnyway -- read file anyway, if np.loadtxt raise ValueError
		READTIMEOUT -- seconds for file reading by hands
		project -- project all particles to sphere with constant radius
		'''
		self.filename=fname
		Norm=1.
		fileheader=""
		if(verbose):print("Reading file",fname,"...")
		with open(fname,"r") as fl: # Norm converts inner parrots to kpc
			tmp=fl.readline()
			import re
			Norm=float(re.compile("T/Mpc = ([0-9\.]+)").findall(tmp)[0])*1e3
			fileheader=tmp
			while("#" in tmp):
				tmp=fl.readline()
				fileheader+=tmp
		self.Norm=Norm
		columns=fileheader.split("\n")[-3].split("\t") # Columns in the file
		self.columns=columns
		start_t=time.time()
		try:
			self.ftable=np.loadtxt(fname) # All file
		except ValueError as e:
			print("Oops... File is broken. NumPy load failed. Try to load what posible")
			if(readAnyway):
				self.ftable=np.loadtxt(fixBitFile(fname))
			else:raise e
			if(verbose):print("Readed",len(self.ftable),"particles")
		if(verbose and time.time()-start_t>15):print("Wait a moment!")
		self.weights=self.ftable[:,1]
		if(verbose):print("Min weight:",self.weights.min(),"Max weight:",self.weights.max())
		# End particle coords
		self.x=self.ftable[:,2]*Norm
		self.y=self.ftable[:,3]*Norm
		self.z=(1-self.ftable[:,4])*Norm
		tmp=(self.x**2+self.y**2+self.z**2)**0.5
		self.maxDistance=max(tmp)
		self.meanDistance=tmp.mean()
		self.minDistance=min(tmp)
		if(verbose):
			print("Max distance (by datafile)",self.maxDistance)
			print("Mean distance (by datafile)",self.meanDistance)
			print("Min distance (by datafile)",self.minDistance)
		if project:
			if type(project)==int or type(project)==float:
				self.meanDistance=project
				print("Mean distance overwrited to",self.meanDistance)
			self.x=self.x/tmp*self.meanDistance
			self.y=self.y/tmp*self.meanDistance
			self.z=self.z/tmp*self.meanDistance
			print("All particles was projected to",self.meanDistance)
		# End particle directions
		self.Px=self.ftable[:,5]
		self.Py=self.ftable[:,6]
		self.Pz=1-self.ftable[:,7]
		self.in_target=np.zeros_like(self.weights).astype(np.bool)
		self.in_target_num=0
		try:
			self.in_target=self.ftable[:,columns.index("in_target")].astype(np.bool)
			self.in_target_num=sum(self.in_target)
			print("All particles:",len(self.in_target),", in target:",self.in_target_num)
		except:
			print("All particles:",len(self.in_target),"No target column found")
		self.params_dict=dict()
		try:
			p=pathlib.Path(fname).resolve().parent.parent.as_posix()
			self.label=p[p.rindex("/")+1:]
			p+="/params"
			params=open(p.replace("//","/"),'r').readlines()
			if(verbose):print("Load run params from",p,"file")
			for l in params:
				if(l.startswith("#")):continue
				if(l.startswith("--")):
					par=l.replace("\n","").split("\t")
					try:
						self.params_dict[par[0][2:]]=par[1]
						if par[0][2:]=="source-x": self.source_x = float(par[1])
						if par[0][2:]=="source-y": self.source_y = float(par[1])
						if par[0][2:]=="source-z": self.source_z = float(par[1])
					except IndexError:
						self.params_dict[par[0][2:]]=""
			self.params_dict['z']=float(self.params_dict['z'])
			if(self.params_dict['z']>1):
				print("Warning: redshift is so big! You cannot use current z->D approximation!")
			redshiftV=3e5*((1+self.params_dict['z'])**2-1)/((1+self.params_dict['z'])**2+1) # km/sec
			H = 70 # km/(sec*Mpc)
			maxDistance=redshiftV/H * 1e3 # in Kpc
			print("Max distance (by redshift):",maxDistance,"Kpc")
			self.paramDistance=maxDistance
		except Exception as e:
			print("Warning: failed to parse run commands:",e)
			if(verbose):raise e
		self.target_x=None
		self.target_y=None
		self.target_z=None
		if("observer-x" in self.params_dict and self.in_target_num>0):
			self.target_x=float(self.params_dict["observer-x"])
			self.target_y=float(self.params_dict["observer-y"])
			self.target_z=float(self.params_dict["observer-z"])
			print("Automatically found target in params file coords:",(self.target_x,self.target_y,self.target_z))
			#self.checkTargetCoords()
			if('source_x' in self.keys()):
				self.calcDirectionOnTarget(verbose=verbose)
			elif(verbose):
				print("Source direction not specified. Restart CRbeam with -os option for more features")
		if(verbose):print("File",fname,"loaded!")
	def __getitem__(self,k):
		'''
		Return param from ../params file if available
		Or return column from datafile
		'''
		try:
			return self.params_dict[k]
		except KeyError:
			return self.ftable[:,self.columns.index(k)]
	def __len__(self):
		'''
		Return number of all particles
		'''
		return len(self.x)
	def keys(self):
		'''
		Return set of available keys
		'''
		return set(self.params_dict)|set(self.columns)
	def calcDirectionOnTarget(self,verbose=True):
		'''
		Calculate (theta,phi) direction for hit in target by --jet-theta and --jet-phi params in CRbeam
		'''
		#from PlotHelpers import xyz2thetaPhi,approxThetaParam
		if(self.in_target.sum()==0):raise Exception("No particles in target")
		tmpx,tmpy,tmpz=self['source_x'][self.in_target],self['source_y'][self.in_target],self['source_z'][self.in_target]
		theta,phi=xyz2thetaPhi(tmpx.mean(),
								tmpy.mean(),
								tmpz.mean())
		if(verbose):
			print("Direction on target: theta =",theta,", phi =",phi)
			rjett=np.pi/2
			if(self.in_target.sum()==1):
				rjett=float(self['observer-radius'])/((tmpx**2+tmpy**2+tmpz**2)**0.5).min()
				rjett=min(2*np.arctan(rjett),np.pi)
				print("Warning: prediction only for one particle!")
			else:
				rjett=approxThetaParam(tmpx,tmpy,tmpz)
			print("> Recommended parameters for next CRbeam run:")
			print("\t-jet",rjett,"-jett",theta,"-jetp",phi)
		self.target_direction=theta,phi
		return theta,phi
	def checkTargetCoords(self):
		'''
		Check target:
		Particles 'in_target' should be consistent with target_x,y,z
		Deprecated.
		'''
		x_in=self.x[self.in_target]
		y_in=self.y[self.in_target]
		z_in=self.z[self.in_target]
		r=max(x_in**2+y_in**2+z_in**2)**0.5
		x0=np.array([x_in.mean(),y_in.mean(),z_in.mean()])
		x0=x0/np.linalg.norm(x0)*r
		if(not self.target_x is None):
			trueRes=np.array([self.target_x,self.target_y,self.target_z]) # Successfully readed from params-file
			if( sum((x0-trueRes)**2)**0.5>r ):
				print("Warning: target coordinats inconsistent. Check it!")
				print("Warning: by file:",[self.target_x,self.target_y,self.target_z])
				print("Warning: by calc:",list(x0))
				return False
		return True
	def manualTarget(self,r,x,y=None,z=None):
		'''
		Reset target coordinates by hands and recalc in_target
		Call: manualTarget(radius, [x,y,z]) or manualTarget(radius,x,y,z)
		'''
		if(y is None and z is None):
			y=x[1]
			z=x[2]
			x=x[0]
		self.target_x,self.target_y,self.target_z=x,y,z
		self.in_target=((self.x-x)**2+(self.y-y)**2+(self.z-z)**2)<r**2
		self.in_target_num=sum(self.in_target)
