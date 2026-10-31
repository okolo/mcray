import FileHelpers
import sys,os
import numpy as np
import re

outputs=re.compile("multiAnother2Void_D[0-9]+_R1000")

dirs=os.listdir("output")

def getNaiveFlux(fl):
 	referentPoint=np.array([float(fl['observer-x']),float(fl['observer-y']),float(fl['observer-z'])])
 	D=sum(referentPoint**2)**0.5
 	return (1-np.cos(np.arctan(float(fl['observer-radius'])/D)))/2

def getFlux(fl):
 	k=2/(1-np.cos(float(fl['jet-angle'])))
 	#print("Factor:",k)
 	N=k*len(fl.in_target)
 	#print("N real:",len(fl.in_target))
 	#print("N isotropic:",N)
 	inp=fl.in_target_num
 	return inp/N

res=open("/tmp/resVoid.dat","w")

print("E","naive","der",file=res)

for dir in [outputs.findall(dir) for dir in dirs]:
	if dir==[]:continue
	dir=dir[0]
	fl=FileHelpers.CRbeamFile("output/"+dir+"/z0/proton")
	naive=getNaiveFlux(fl)
	realfl=getFlux(fl)
	der=realfl/naive
	D=re.compile("multiAnother2Void_D([0-9]+)_R1000").findall(dir)[0]
	print(D,naive,der,file=res)

res.close()
