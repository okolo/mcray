#!/usr/bin/python3
#!/usr/bin/env python3
# Requires:
# sudo apt install python3-pip python3-numpy python3-matplotlib python3-tk
# /usr/bin/python3 -m pip install healpy
# /usr/bin/python3 -m pip install tqdm
# Anaconda is not recommended
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
import tkinter as tk
from tkinter import messagebox

plt.ion()

ENABLE_FLIP=False#True # направить ось Z на нас, DEPRECATED: теперь используется rot, UPD: в healpy баг, лучше делать поворот своими руками.
rot=(0,90,0) # Поворот картинки осью z к центру, NOTE: в healpy баг, осторожно.
enableLog=True # Логарифмический масштаб картинки, UPD: используются штатный параметр healpy
enableRef=True # Отметить центр/наблюдателя точкой
enableGrad=True # Рисовать сетку
#enableBigTargetFix=False # Исправлять размазывание источника с точки зрения большого наблюдателя. Не влияет, если наблюдатель далеко
nside=256

fname="photon"
try:
	fname=sys.argv[1]
except:
	print("No file specified, try to open \"photon\"")

print("Loading file...")
Norm=1.
fileheader=""
with open(fname,"r") as fl: # Коэффициент преобразует внутренние попугаи в килопарсеки
	tmp=fl.readline()
	import re
	Norm=float(re.compile("T/Mpc = ([0-9\.]+)").findall(tmp)[0])*1e3
	fileheader=tmp
	while("#" in tmp):
		tmp=fl.readline()
		fileheader+=tmp

dirname = fname[:fname.rindex("/z0/")]

def getParamByName(dirname,paramname):
	params=open(dirname+"/params","r").readlines()
	for line in params:
		if line.startswith(paramname):
			try:
				return float(line.split("\t")[1])
			except:
				return line.split("\t")[1]
	return None

sourceX = getParamByName(dirname,"--source-x")
sourceY = getParamByName(dirname,"--source-y")
sourceZ = getParamByName(dirname,"--source-z")

columns=fileheader.split("\n")[-3].split("\t")

import time
start_t=time.time()
phot=np.loadtxt(fname)
if(time.time()-start_t>15):print("Wait a moment!")
x=phot[:,2]*Norm# - getParamByName(dirname,"--source-x")
y=phot[:,3]*Norm# - getParamByName(dirname,"--source-y")
z=(1-phot[:,4])*Norm# - getParamByName(dirname,"--source-z")

try:
	x-=getParamByName(dirname,"--source-x")
	y-=getParamByName(dirname,"--source-y")
	z-=getParamByName(dirname,"--source-z")
except Exception as e:
	print(e)

tmp=(x*x+y*y+z*z)**0.5
print("Max distance (by datafile)",max(tmp))
print("Min distance (by datafile)",min(tmp))
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
	in_target=phot[:,columns.index("in_target")].astype(np.bool)
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
	target_x=float(params_dict["observer-x"])-float(params_dict["source-x"])
	target_y=float(params_dict["observer-y"])-float(params_dict["source-y"])
	target_z=float(params_dict["observer-z"])-float(params_dict["source-z"])
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
	if(lbl==""):lbl=sys.argv[1]
except:
	lbl=sys.argv[1]

def xyz2thetaPhi(x,y,z):
	'''
	Преобразует переданные декартовы координаты в углы на сфере в радианах
	Принимает числа или NumPy-одномерные массивы одинаковой длины
	Меняет X и Z местами, если ENABLE_FLIP истина
	'''
	if(ENABLE_FLIP):
		x,z=z,x
		z*=-1 # Чтобы совпасть с багнутым вариантом, который кроме бага идеологически верен
	r=(x*x+y*y+z*z)**0.5
	theta=np.arccos(z/r)
	phi=1.5707963267948957
	try:
		phi=np.arctan(y/x)+(x<0)*3.1415926535897932
		phi+=(phi<0)*2.*3.1415926535897932
	except:
		print("Warning: singular phi")
	return theta,phi

def thetaPhi2xyz(theta,phi,R):
	'''
	Transform spherical coords in radians to x,y,z coordinates
	Input: theta, phi, radius
	Return x,y,z
	'''
	Tx=R*np.sin(theta)*np.cos(phi)
	Ty=R*np.sin(theta)*np.sin(phi)
	Tz=R*np.cos(theta)
	return Tx,Ty,Tz

def calcTargetCoords(x_in,y_in,z_in,verbose=True):
	'''
	Возвращает координаты центра сферы по заданным массивам трёхмерных координат.
	Если координаты наблюдателя были получены из парам-файла, то возвращает их без вычислений
	'''
	global target_x,target_y,target_z
	r=max(x_in*x_in+y_in*y_in+z_in*z_in)**0.5
	x0=np.array([x_in.mean(),y_in.mean(),z_in.mean()])
	x0=x0/np.linalg.norm(x)*r
	if(not target_x is None):
		trueRes=np.array([target_x,target_y,target_z]) # Успешно стырено из парам-файла
		if( sum((x0-trueRes)**2)**0.5>r ):
			print("Warning: target coordinats inconsistent. Check params file!")
			print("Warning: by file:",[target_x,target_y,target_z],"(it will use later)")
			print("Warning: by calc:",list(x0))
		return trueRes
	return x0 # Похоже, это самый нормальный вариант угадывания, не внося кардинальных изменений в формат вывода файлов...
	# в центре сферы расстояния до всех точек равны + запрет отходить дальше чем на разброс по координатам
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
	source_x=phot[:,columns.index("source_x")]*Norm
	source_y=phot[:,columns.index("source_y")]*Norm
	source_z=phot[:,columns.index("source_z")]*Norm
	theta,phi=xyz2thetaPhi(source_x.mean(),source_y.mean(),source_z.mean())
	print("Jet direction: theta =",theta,", phi =",phi)
	#theta,phi=xyz2thetaPhi(source_x[in_target].mean(),source_y[in_target].mean(),source_z[in_target].mean())
except:
	print("Source direction are not specified")


import healpy as hp
from tqdm import tqdm

thetas,phis=xyz2thetaPhi(x,y,z)
Dmean=((x**2+y**2+z**2)**0.5).mean()
print("Force reset all D to ",Dmean)
x,y,z=thetaPhi2xyz(thetas,phis,Dmean)


def hplot(thetas,phis,weights,postprocess=lambda x:x,enableRef=enableRef):
	'''
	Plot thetas,phis angle array by HealPy
	'''
	npix = hp.nside2npix(nside)
	# Go from HEALPix coordinates to indices
	indices = hp.ang2pix(nside, thetas, phis)
	# Initate the map and fill it with the values
	hpxmap = np.zeros(npix, dtype=np.float)
	tmp=weights
	for i in tqdm(range(len(tmp))):
		hpxmap[indices[i]] += tmp[i]
	hpxmap=postprocess(hpxmap)
	if(enableLog):
		#hpxmap=np.log(1.+hpxmap)
		hpxmap+=1.
		norm='log'
	else:
		norm=None
	print("DEBUG 216: max(hpxmap)=",max(hpxmap))
	# Inspect the map
	xsize=max(round(1.05*(len(hpxmap)*2/3.14159)**0.5)*2,800)
	print("max(hpxmap) =",max(hpxmap))
	hp.mollview(hpxmap,rot=rot,norm='log',xsize=xsize)
	savepath=fname[:fname.rindex("/z0/")]
	np.savetxt(savepath+"/general_view.hpmap",hpxmap)
	if(enableGrad):hp.graticule()
	if(enableRef):
		t,p=xyz2thetaPhi(0., 0., 1.)
		if(sum(in_target)>0):
			target_coords=calcTargetCoords(x[in_target],y[in_target],z[in_target],verbose=False)
			t,p=xyz2thetaPhi(target_coords[0],target_coords[1],target_coords[2])
			if(np.isnan(p) or np.isinf(p)):p=0.0
			if(np.isnan(t) or np.isinf(t)):t=0.0
		print("Mark target(if exist) as red point (theta,phi) =",(t,p))
		hp.projscatter(t, p,s=250,color='red',marker='.')
	return hpxmap

hplot(thetas,phis,weights,enableRef=sum(in_target)>0)
plt.title("End particles sphere by source "+lbl)
plt.savefig(fname="/tmp/MyFigure"+str(nside)+".png",dpi=240,orientation='landscape')
plt.show(block=False)

def showTargetFeatures(target_coords,enablePlots=[False,True,True],referentPoint=None):
	'''
	Show 3 plots by np.array([x,y,z]) coordinates of target
	in_target array should be not trivial and has the same dimensions as x,y,z and etc
	'''
	if(enablePlots[0]):
		'''
		Куда надо помещать детекторы, чтобы увидеть частицы?
		Бессмысленный график, отключён по умолчанию
		'''
		thetas,phis=xyz2thetaPhi(target_coords[0]-x[in_target],target_coords[1]-y[in_target],target_coords[2]-z[in_target])
		hplot(thetas,phis,weights[in_target],enableRef=False)
		plt.title("Where on the target observer see "+lbl)
		plt.show(block=False) ######### PLOT
	if(enablePlots[1]):
		'''
		Куда источник должен пускать частицы, чтобы попасть в цель?
		Если при запуске CRbeam не был указан спец.ключ, то вываливает ошибку
		'''
		try:
			# Source directions
			thetas,phis=xyz2thetaPhi(source_x,source_y,source_z)
			hplot(thetas,phis, 4000.*in_target-0.3, postprocess=lambda x:(x!=0).astype(np.float)+2.*(x>0).astype(np.float) )
			plt.title("Initial particles directions "+lbl)
			plt.show(block=False) ######### PLOT
		except Exception as e:
			print(e,"May be source direction not specified?")
	if(enablePlots[2]):
		'''
		Что наблюдатель видит на небе?
		'''
		# Direction
		if(referentPoint is None):
			r=(x[in_target]**2+y[in_target]**2+z[in_target]**2)**0.5
	#		referentPoint=np.array([0,0,1])*r.max()
			thetas,phis=xyz2thetaPhi(Px[in_target]-x[in_target]/r,Py[in_target]-y[in_target]/r,Pz[in_target]-z[in_target]/r)
			enableBigTargetFix = False
			if(enableBigTargetFix): # DEPRECATED
				print("Warning: enableBigTargetFix is bad deprecated param!")
				print("Fix directions...")
				thetasFix,phisFix=xyz2thetaPhi(x[in_target],y[in_target],z[in_target])
				thetas-=thetasFix
				phis-=3.1415926535897932*(thetas<0) # Перескок по тетте -- поворот на 180 градусов
				phis+=2*3.1415926535897932*(phis<0)
				thetas=np.abs(thetas)
#				#phis-=phisFix
#				#print("Fixed thetas: ",thetas.min(),thetas.max())
		if(referentPoint is not None):
			'''
			Нетривиальная математика:
			1)Берём направление на источник e1=r=(x,y,z)
			2)Берём направление на referentPoint: e2=referentPoint-r;
				e2=e2-r*(r,e2)/(r,r) -- строим ортогональный базис
			3)Правая тройка векторов, в студию: e3=[e1,e2], вращайте базис!
			4)Новые координаты для построения:
				p_x=(P,e1)
				p_y=(P,e2)
				p_z=(P,e3)
			'''
			r=np.vstack((x[in_target],y[in_target],z[in_target]))
			r0 = np.zeros_like(r)
			if source_x is not None:
				r0 = np.vstack((source_x[in_target],source_y[in_target],source_z[in_target]))
				r = r-r0
			e1=r/np.linalg.norm(r,axis=0) # This minus by hands for observer in (0,0,0) and random source
			e2=np.vstack((referentPoint[0]-(x[in_target]-source_x[in_target]),
							referentPoint[1]-(y[in_target]-source_y[in_target]),
							referentPoint[2]-(z[in_target]-source_z[in_target])))
			e2-=e1*(e1*e2).sum(axis=0)
			e2/=np.linalg.norm(e2,axis=0)
			e3=np.cross(e1,e2,axisa=0,axisb=0,axisc=0)
			P=np.vstack((Px[in_target],Py[in_target],Pz[in_target]))
			# Проекция на е1 ставится как Z, чтобы направление на источник было сонаправлено с осью Z
			thetas,phis=xyz2thetaPhi((e3*P).sum(axis=0),(e2*P).sum(axis=0),(e1*P).sum(axis=0))
		hpxmap=hplot(thetas,phis,weights[in_target],enableRef=False)
		savepath=fname[:fname.rindex("/z0/")]
		print("savepath:",savepath)
		np.savetxt(savepath+"/target_view.hpmap",hpxmap)
		# TODO: Дописать максимальный угол
		theta_max=thetas.argmax()
		phi_max=phis[theta_max]+0.01
		theta_max=thetas[theta_max]
		theta_max=int(round(theta_max/3.14159*180))
		theta_max=theta_max*3.14159/180
		max_angle_color='orange'
		for i in range(64):
			hp.projscatter(theta_max*1.005, phi_max+i/64*6.28,s=15,color=max_angle_color,marker='.')
		hp.visufunc.projtext(theta_max*1.02, phi_max,
							color=max_angle_color,
							s='$\\theta='+str(int(round(theta_max/3.14159*180)))+'^\\circ$',
							fontsize=24)
		if(referentPoint is not None):
			t,p=xyz2thetaPhi(0., 1., 0.)
			hp.projscatter(t, p,s=70,color='orange',marker='.')
		plt.title("")#"What observer see "+lbl)
	plt.show(block=False)

def approxThetaParam(x,y,z):
	'''
	Находит оценку сверху для угла theta ширины джета при симуляциях, чтобы попасть в наблюдателя
	Принимает NumPy-одномерные массивы x,y,z одинаковой длины
	'''
	# TODO: Придумать способ учёта частиц, которые НЕ попали в цель: посчитать угол раствора до ближайшей не попавшей частицы
	r=(x*x+y*y+z*z)**0.5
	R=r.min()
	dx=x.max()-x.min()
	dy=y.max()-y.min()
	dz=z.max()-z.min()
	dr=max(dx,dy,dz)
	return np.arctan(dr/R)

def getNaiveFlux(D,R):
	return (1-np.cos(np.arctan(R/D)))/2

def getFlux(jet_angle,in_target,particles):
	k=2/(1-np.cos(jet_angle))
	#print("Factor:",k)
	N=k*particles
	#print("N real:",len(fl.in_target))
	#print("N isotropic:",N)
	return jet_angle/N

if(sum(in_target)>0):
	print("\n================ Found target! ================")
	target_coords=calcTargetCoords(x[in_target],y[in_target],z[in_target])
	print("> Target coordinates:",target_coords)
	print("> Catched:",sum(in_target),"("+str(round(sum(in_target)*100/len(in_target),5))+"%)")
	if(sum(in_target)>2 and source_x is not None):
		rjett=approxThetaParam(source_x[in_target],source_y[in_target],source_z[in_target])
		rtheta,rphi=xyz2thetaPhi(source_x[in_target].mean(),source_y[in_target].mean(),source_z[in_target].mean())
		print("> Recommended parameters for next CRbeam run:")
		print(f"\t-jet {rjett} -jett {rtheta} -jetp {rphi}")
	print("===============================================")
	showTargetFeatures(target_coords)

def latlong2xyz(lat,long,R=None):
	'''
	Transform lat&long coords in degrees to x,y,z coordinates
	Input: latitude, longitude, radius (optional)
	Return np.array([x,y,z])
	'''
	phi=long/180*3.1415926
	if(phi<0):phi+=2*3.1415926
	theta=(90-lat)/180*3.1415926
	if(R is None):
		R=max((x*x+y*y+z*z)**0.5)
	Tx=R*np.sin(theta)*np.cos(phi)
	Ty=R*np.sin(theta)*np.sin(phi)
	Tz=R*np.cos(theta)
	return np.array([Tx,Ty,Tz])

class Application(tk.Frame):
	def __init__(self, master=None,r=100):
		super().__init__(master)
		self.master = master
		self.master.title("CRbeam Plotter")
		if("anaconda" in tk.__file__.lower()): # font bug
			self.master.minsize(950, 393)
		else:
			self.master.minsize(623, 393)
		self.pack(fill="both")
		self['bg']="#f9f9f9"
		self.create_widgets(r)
		self.blocked=False
	def create_widgets(self,r):
		# First line
		currow=1
		tk.Label(self,text="Параметры наблюдателя:",font='Arial 14',bg=self['bg']).grid(row = currow, column = 1,columnspan = 2)
		# Second line
		currow=2
		tk.Label(self,text="       lat=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 1,columnspan = 1)
		tk.Label(self,text="deg"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 3,columnspan = 1)
		self.mlat = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.mlat.insert(0,"90.0")
		self.mlat.grid(row = currow, column = 2,columnspan = 1)
		# Therd line
		currow=3
		tk.Label(self,text="    long=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 1,columnspan = 1)
		tk.Label(self,text="deg"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 3,columnspan = 1)
		self.mlong = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.mlong.insert(0,"0.0")
		self.mlong.grid(row = currow, column = 2,columnspan = 1)
		# 4-th line
		currow=4
		tk.Label(self,text="D(obs)=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 1,columnspan = 1)
		tk.Label(self,text="kpc"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 3,columnspan = 1)
		self.mD = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.mD.insert(0,str(round(max((x*x+y*y+z*z)**0.5),1)))
		self.mD.grid(row = currow, column = 2,columnspan = 1)
		# 5-th line
		currow=5
		tk.Label(self,text="R(obs)=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 1,columnspan = 1)
		tk.Label(self,text="kpc"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 3,columnspan = 1)
		self.mR = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.mR.insert(0,str(r))
		self.mR.grid(row = currow, column = 2,columnspan = 1)
		# 6-th line
		currow=6
		tk.Label(self,text="Параметры Healpy:",font='Arial 13',bg=self['bg']).grid(row = currow, column = 1,columnspan = 2)
		# 7-th line
		currow=7
		tk.Label(self,text="nside=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 1,columnspan = 1)
		self.mNside = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.mNside.insert(0,str(nside))
		self.mNside.grid(row = currow, column = 2,columnspan = 1)
		#self.mNside.bind("<Enter>", self.setObs)
		# 8-th line
		currow=8
		tk.Label(self,text="Строить графики:",font='Arial 13',bg=self['bg']).grid(row = currow, column = 1,columnspan = 2)
		# 9-th line
		currow=9
		self.firstPlot=tk.BooleanVar(self)
		tk.Checkbutton(self, text="Распределение частиц вокруг источника    ",font='Arial 12',bg=self['bg'],onvalue=1,offvalue=0,variable=self.firstPlot).grid(row = currow, column = 1,columnspan = 3)
		self.firstPlot.set(1)
		# 10-th line
		currow=10
		self.secondPlot=tk.BooleanVar(self)
		tk.Checkbutton(self, text="Куда попадают частицы на наблюдателе   ",font='Arial 12',bg=self['bg'],onvalue=1,offvalue=0,variable=self.secondPlot).grid(row = currow, column = 1,columnspan = 3)
		self.secondPlot.set(0)
		# 11-th line
		currow=11
		self.therdPlot=tk.BooleanVar(self)
		self.therdPlot.set(0)
		if(source_x is None):
			pass
		tk.Checkbutton(self, text="Куда надо пустить частицы для попадания",font='Arial 12',bg=self['bg'],onvalue=1,offvalue=0,variable=self.therdPlot).grid(row = currow, column = 1,columnspan = 3)
		# 12-th line
		currow=12
		self.fourthPlot=tk.BooleanVar(self)
		self.fourthPlot.set(1)
		tk.Checkbutton(self, text="Что видит наблюдатель на небе"+16*" ",font='Arial 12',bg=self['bg'],onvalue=1,offvalue=0,variable=self.fourthPlot).grid(row = currow, column = 1,columnspan = 3)
		# 13-th line
		currow=13
		self.clearPlot=tk.BooleanVar(self)
		self.clearPlot.set(1)
		tk.Checkbutton(self, text="Очищать все рисунки перед построением"+2*" ",font='Arial 12',bg=self['bg'],onvalue=1,offvalue=0,variable=self.clearPlot).grid(row = currow, column = 1,columnspan = 3)
		# 14-th
		currow=14
		self.setObsBtn = tk.Button(self,text="Построить",font='Arial 12',command=self.setObs)
		self.setObsBtn.grid(row = currow, column = 1,columnspan = 1)
		self.stateL=tk.Label(self,text="Готов",font='Arial 12',bg=self['bg'])
		self.stateL.grid(row = currow, column = 2,columnspan = 3)
		# Output info
		currow=15
		self.outInfo=tk.Label(self,text="",font='Arial 12',bg=self['bg'])
		self.outInfo.grid(row = currow, column = 1,columnspan = 3)
		# Reference point
		currow=1
		tk.Label(self,text="Референтная точка:",font='Arial 14',bg=self['bg']).grid(row = currow, column = 4,columnspan = 2)
		# Second line
		currow=2
		tk.Label(self,text="       lat=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 4,columnspan = 1)
		tk.Label(self,text="deg"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 6,columnspan = 1)
		self.reflat = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.reflat.insert(0,"90.0")
		self.reflat.grid(row = currow, column = 5,columnspan = 1)
		# Therd line
		currow=3
		tk.Label(self,text="    long=",font='Arial 12',bg=self['bg']).grid(row = currow, column = 4,columnspan = 1)
		tk.Label(self,text="deg"+6*" ",font='Arial 12',bg=self['bg']).grid(row = currow, column = 6,columnspan = 1)
		self.reflong = tk.Entry(self,font='Arial 12',bg=self['bg'],width=10)
		self.reflong.insert(0,"0.0")
		self.reflong.grid(row = currow, column = 5,columnspan = 1)
	def setObs(self,param=None):
		global nside,target_coords,in_target
		if(self.blocked):return
		self.blocked=True # Блокирует попытки нажать кнопку ещё раз
		try:
			self.stateL['text']='Считается...'
			self.update()
			lat=float(self.mlat.get().replace(",","."))
			long=float(self.mlong.get().replace(",","."))
			r=float(self.mR.get().replace(",","."))
			D=float(self.mD.get().replace(",","."))
			flux_avg=getNaiveFlux(D,r)
			print("Flux average:",flux_avg)
			nside=int(self.mNside.get().replace(",","."))
			print("Set observer with params:",lat,long,r,nside)
			needPlots=[self.firstPlot.get(),self.secondPlot.get(),self.therdPlot.get(),self.fourthPlot.get()]
			target_coords=latlong2xyz(lat,long,D)
			target_x=target_coords[0] # Dirty trick!
			target_y=target_coords[1]
			target_z=target_coords[2]
			self.outInfo['text']='Расстояние до наблюдателя:\nD='+str(round(sum(target_coords**2)**0.5,1))+' kpc\n'
			self.outInfo['text']+='\nКоординаты наблюдателя (kpc):\nx='+str(round(target_x+sourceX,2))+' y='+str(round(target_y+sourceY,2))+' z='+str(round(target_z+sourceZ,2))
			self.update()
			in_target=((target_coords[0]-x)**2+(target_coords[1]-y)**2+(target_coords[2]-z)**2)<r**2
			self.outInfo['text']+='\n\nПопало в цель '+str(sum(in_target))+' частиц'
			if(self.clearPlot.get()):
				plt.close('all')
			self.update()
			reflat=float(self.reflat.get().replace(",","."))
			reflong=float(self.reflong.get().replace(",","."))
			refPoint=latlong2xyz(reflat,reflong,D)
			self.stateL['text']='Строится...'
			self.update()
			if(needPlots[0]):
				hplot(thetas,phis,weights)
				if(refPoint is not None):
					t,p=0,0
					if(reflat!=90 or reflong!=0):
						t,p=xyz2thetaPhi(refPoint[0], refPoint[1], refPoint[2])
					hp.projscatter(t, p,s=70,color='orange',marker='.')
				plt.title("")#"End particles sphere by source "+lbl)
				plt.show(block=False)
			if(sum(in_target)>0):
				print("Real flux:",sum(in_target)/len(in_target))
				print("F/F_{av}",sum(in_target)/len(in_target)/flux_avg)
				print("===============================================")
				print("> Target coordinates:",target_coords)
				print("> Catched:",sum(in_target),"("+str(round(sum(in_target)*100/len(in_target),5))+"%)")
				recomD=(x[in_target]**2+y[in_target]**2+z[in_target]**2)**0.5
				recomD=recomD.mean()
				print("Recomended D for GUI:",recomD)
				if(sum(in_target)>2 and source_x is not None):
					rjett=approxThetaParam(source_x[in_target],source_y[in_target],source_z[in_target])
					rtheta,rphi=xyz2thetaPhi(source_x[in_target].mean(),source_y[in_target].mean(),source_z[in_target].mean())
					print("> Recommended parameters for next CRbeam run:")
					print(f"\t-jet {rjett} -jett {rtheta} -jetp {rphi}")
					print(f"\t-obs 1 -obsx {target_x+sourceX} -obsy {target_y+sourceY} -obsz {target_z+sourceZ}")
				print("===============================================")
				showTargetFeatures(target_coords,needPlots[1:],referentPoint=refPoint)
			else:
				print("При заданных параметрах в цель ничего не попало")
			self.stateL['text']="Готово"
			self.update()
		except Exception as e:
			self.stateL['text']="Ошибка: "+str(e)
			self.update()
			raise e
		self.blocked=False

def GUIloop():
	root = tk.Tk()
	app = Application(master=root)
	app.mainloop()

GUIloop()
