#!/net/minus/home/dolgikh/anaconda3/bin/python3

from FileHelpers import *
import psutil,random
import os,time,sys

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

def getThetaByR(R,Robs=200):
	'''
	Convert distance between source and radius to approxjet angle in degrees
	R - distance (kpc)
	Robs - observer's radius (kpc)
	'''
	model1 = 0.1336*R**0.49
	model2 = 0.0882*R**0.53
	if isinstance(R,(float,int)):
		return max(model1,model2)*1.55+(Robs/(R+0.01))/3.14159*180
	sel = (model2>model1).astype(float)
	return (sel*model2+(1-sel)*model1)*1.55+(Robs/(R+0.01))/3.14159*180

# The center
obsx = 0
obsy = 0
obsz = 0

obsr = 200 # kpc
D = 100*1000 # kpc, max distance
Lcorr = 1 # Mpc

'''def getRandCoords():
	x,y,z = D*(1-2*random.random()),D*(1-2*random.random()),D*(1-2*random.random())
	while x**2 + y**2 + z**2>D**2:
		x,y,z = D*(1-2*random.random()),D*(1-2*random.random()),D*(1-2*random.random())
	return x,y,z'''

def getRandCoords():
	x,y,z = D*(1-2*random.random()),D*(1-2*random.random()),D*(1-2*random.random())
	while x**2 + y**2 + z**2>D**2 or x**2 + y**2 + z**2<(50*1000)**2: # Мы с Басей настолько устали, что я не хочу думать над чем-то более умным
		x,y,z = D*(1-2*random.random()),D*(1-2*random.random()),D*(1-2*random.random())
	return x,y,z

N=11*100*1000

def runCRbeam():
	z=round(D2z(obsr+(source_x**2+source_y**2+source_z**2)**0.5 )*1.3,7)
	nth=0
	for i in range(5):
		time.sleep(2)
		nth+=int((1-psutil.cpu_percent()/100)*psutil.cpu_count())-1
	nth=int(nth/5)
	if(nth<1):nth = 1
	if(nth>60):print("Cluster is free, very nice!")
	if(nth>110):nth = 110
	cmd="./CRbeam -p 10 -z "+str(z)+" -E 1e19 -e 1e15 -mft -mf 1e-9 -nth "+str(nth)+" -os -o output/"
	cmd+="srcRand42L_"+str(source_x)+"_"+str(source_y)+"_"+str(source_z)+"_R"+str(obsr)
	cmd+=" -obs 1 -obsr "+str(obsr)
	cmd+=" -jet "+str(jet_angle)+" -jett "+str(jet_theta)+" -jetp "+str(jet_phi)
	cmd+=" -N "+str(N)+" -rseed 42 --lEGMF 1"
	cmd+=" -srcx "+str(source_x)+" -srcy "+str(source_y)+" -srcz "+str(source_z)
	print("RUN COMMAND: "+cmd)
	os.system(cmd)

while True:
	try:
		source_x,source_y,source_z = getRandCoords()
		source_x = round(source_x,1)
		source_y = round(source_y,1)
		source_z = round(source_z,1)
		jet_theta,jet_phi=xyz2thetaPhi(-source_x,-source_y,-source_z)
		jet_angle=getThetaByR( (source_x**2+source_y**2+source_z**2)**0.5, obsr )/180*3.1416
		runCRbeam()
		# sys.exit(0)
	except KeyboardInterrupt:
		print("Выходим")
		sys.exit(0)
