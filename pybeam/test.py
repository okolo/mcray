#!/usr/bin/env python3
import pybeam
import numpy as np

pybeam.testCPP()

# pybeam.testMCpropag()

pybeam.cosmologyInit(0.02)

print("pybeam.units.eV =",pybeam.units.eV)
print("pybeam.units.Mpc =",pybeam.units.Mpc)

cmbTemp = 2.73/pybeam.units.phTemperature_mult/pybeam.units.Eunit

print("cmbTemp =",cmbTemp)

pOutput = pybeam.RawOutput3D("/tmp/test_out",True,False,False,False,False,False)
pOutput.SetOutputDir("/tmp/test_out/z0",True)
print("Created RawOutput3D:",pOutput)

backgr = pybeam.CompoundBackground()

print("Created CompoundBackground:",backgr)

b1 = pybeam.PlankBackground(cmbTemp,1e-3*cmbTemp, 1e3*cmbTemp, 0, 1.02)
print("Created PlankBackground:",b1)

backgr.AddComponent(b1)
print("Updated CompoundBackground:",backgr)
backgr.AddComponent(pybeam.MatrixBackground("Inoue12Baseline", "EBL_Inoue/EBL_proper_baseline.dat", False, True),0.5)
print("Updated CompoundBackground:",backgr)

backgr2 = pybeam.CuttedBackground(backgr, 0, 100)
print("Cutted CompoundBackground:",backgr2)

dens = pybeam.BackgroundUtils.CalcIntegralDensity(b1)*pybeam.units.cm3
print("Planck integral density [cm^{-3}]:", dens)

dens = pybeam.BackgroundUtils.CalcIntegralDensity(backgr)*pybeam.units.cm3
print("EBL integral density [cm^{-3}]:", dens)

pybeam.BackgroundUtils.Print(backgr,0,"/tmp/test_out/backgr",100)
print("EBL written to /tmp/test_out/backgr")

logStepK = 10**0.05
fZmax = 0.02
stepZ = fZmax/2
epsRel = 1e-3
centralE1 = 6.3e-10/pybeam.units.Eunit
n1 = 413*pybeam.units.Vunit
centralE1 = 6.3e-10/pybeam.units.Eunit

k1 = pybeam.ConstFunction(centralE1)
print("k1 is ConstFunction:",k1)
print("k1(123) =",k1(123),", centralE1 =",centralE1,", float(k1) =",float(k1))

c1 = pybeam.ConstFunction(n1);
backgrI = pybeam.MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel)
print("MonochromaticBackgroundIntegral from ConstFunctions:",backgrI)

backgrI = pybeam.MonochromaticBackgroundIntegral(n1, centralE1, logStepK, epsRel)
print("MonochromaticBackgroundIntegral from doubles:",backgrI)

backgrI = pybeam.ContinuousBackgroundIntegral(backgr, stepZ, logStepK, fZmax, epsRel)
print("ContinuousBackgroundIntegral from doubles:",backgrI)

rand = pybeam.Randomizer(42)
print("Random seed=42 tests:",rand.Rand(),rand.Rand(),rand.Rand())
print("Saved copy from prev:",0.7275958061218262, 0.8862816691398621, 0.9642069339752197)
print("rand.CreateIndependent() ->",rand.CreateIndependent())

rand = pybeam.Randomizer()
print("Random seed from timestamp tests:",rand.Rand(),rand.Rand(),rand.Rand())

import numpy as np
resultFilter = pybeam.EnergyBasedFilter(10**10*pybeam.units.eV,np.pi)
print("EnergyBasedFilter:",resultFilter)
result = pybeam.Result(resultFilter,True)

result.SetAutoFlushInterval(10)
result.EnableFlushOnCtrlC()
result.AddOutput(pOutput)

print("Result:",result)

particles = pybeam.ParticleStack()
print("ParticleStack:",particles)

rand = pybeam.Randomizer(42)

pe = pybeam.PropagationEngine(particles, result, rand.CreateIndependent())
print("Engine pe:",pe)

alphaThinning = 0 # alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
EThinning = 0*pybeam.units.eV;
thinning = pybeam.PhotonThinning(alphaThinning, EThinning)
print("PhotonThinning:",thinning)
pe.SetThinning(thinning)

discardedParticles = pybeam.ParticleTypeFilter([0,1,2]) # Electron, Positron, Photon
print("discardedParticles:",discardedParticles)

pe.SetPropagationFilter(discardedParticles)

ics = pybeam.ICSInteraction(backgrI, EThinning)
gzk = pybeam.GZK(backgrI)
print("ICSInteraction:",ics)
print("GZK:",gzk)

pe.AddInteraction(ics)
pe.AddInteraction(gzk)
pe.AddInteraction(pybeam.NeutronDecay())
pe.AddInteraction(pybeam.ProtonPPcel(backgrI))

mf = pybeam.TurbulentMF(rand, 1, 1e-9, 1.01, 0.1)
print("TurbulentMF:",mf)
print("TurbulentMF -> GetValueGauss([0,0,0],0)",mf.GetValueGauss([0,0,0],0))
B=mf.GetValueGauss(np.array([-9.00998,-9.00998,-9.00998])*pybeam.units.Mpc,0)
Bref=[-3.71297e-10, 2.47636e-10, 7.13073e-10]
print("TurbulentMF -> GetValueGauss([-9.00998,-9.00998,-9.00998],0)", B )
print("Original value from cluster file", Bref )

diff=(B[0]-Bref[0])**2+(B[1]-Bref[1])**2+(B[2]-Bref[2])**2
diff/=(Bref[0])**2+(Bref[1])**2+(Bref[2])**2
print("Difference: ",diff**0.5)
if diff>1e-4:
	print("TurbulentMF TEST FAILED!")
	exit(1)

print("Save mf...")
mf.savetxt("/tmp/magneticField.out",5,0.5)
print("Saved! You can check it here: /tmp/magneticField.out")

mf2 = pybeam.SavedMF("/tmp/magneticField.out")
Bsaved = mf2.GetValueGauss(np.array([-5,-5,-4.5])*pybeam.units.Mpc,0)
Bturb = mf.GetValueGauss(np.array([-5,-5,-4.5])*pybeam.units.Mpc,0)
diff=(Bsaved[0]-Bturb[0])**2+(Bsaved[1]-Bturb[1])**2+(Bsaved[2]-Bturb[2])**2
diff/=(Bturb[0])**2+(Bturb[1])**2+(Bturb[2])**2
print("Saved mf test:",diff)
if diff>1e-10:
	print("SavedMF TEST FAILED!")
	exit(1)

defl = pybeam.Deflection3D(mf,10)
pe.AddInteraction(defl)

t = pybeam.CosmoTime(0.0085871)
d = -t.t()/pybeam.units.kpc
if abs(d-36643)/36643>0.01:
	print("CosmoTime TEST FAILED!")
	exit(1)

t0 = pybeam.CosmoTime() # Current time
if t0<t:
	print("CosmoTime TEST2 FAILED!")
	exit(1)

D = -pybeam.CosmoTime(0.0024).t()/pybeam.units.kpc # D = 10260.737417323564
print("D =",D)

tEnd = pybeam.CosmoTime()
tEnd.setZ(0)

result.SetEndTime(tEnd)

tStart = pybeam.CosmoTime()
tStart.setZ(fZmax)

dt = tEnd.t()-tStart.t()
t0 = tStart.t()

fBatchSize = 40

particle = pybeam.Particle(pybeam.particle_types.Proton, 0.0024)
print("particle =",particle)

p2 = particle.copy()

if id(particle)==id(p2):
	print("PARTICLE COPY TEST FAILED!")
	exit(1)

particle.Energy = 1e19*pybeam.units.eV

if particle.Energy < 1:
	print("PARTICLE CHANGE PROPERTY TEST FAILED!")
	exit(1)

particle.Time = pybeam.CosmoTime(0.0326962)

if particle.Time == p2.Time:
	print("PARTICLE COPY TEST2 (with time) FAILED!")
	exit(1)

print("particle =",particle)

x = particle.getX()
print("particle X =",x)

particle.setX(np.array([0,0,1])*pybeam.units.kpc*D)

x = particle.getX()
print("Update: particle.getX() = ",x)
print("[0,0,D]*kpc=",np.array([0,0,1])*pybeam.units.kpc*D)

particle.PropagateFreely(pybeam.CosmoTime(0.0024).t())

x = particle.getX()
print("PropagateFreely: particle X =",x, "(must be zeros)")
if (np.array(x)**2).sum()>1e-12:
	print("PARTICLE PropagateFreely TEST FAILED!")
	exit(1)

print("getP ->", particle.getP())
print("getStartP ->", particle.getStartP())
P = np.array([-0.971409,-0.0895426,0.219881])
particle.setP(P,True)
print("getP ->", particle.getP())
print("getStartP ->", particle.getStartP())
print("Input param:", P)

particle.setX(np.array([96648.3,-5035.9,-5601.2])*pybeam.units.kpc)
particle.Time = pybeam.CosmoTime(0.04)
particle.setTracer("/tmp/trajectory.dat")
print("Write trajectory to /tmp/trajectory.dat")

dx = mf.MinVariabilityScale(particle.Time)/10
deltaT = -pybeam.CosmoTime(0.04).t()
print("Time Z ->", particle.Time.z())
print("Deflecting ... For about",round(deltaT/dx + 1.),", kill me hard if it is negative")

import time
start_t = time.time()
defl.Propagate(deltaT,particle,rand)
end_t = time.time()
print("Deflection time:",end_t-start_t)

print("getP ->", particle.getP())
print("getStartP ->", particle.getStartP())
print("getX ->", np.array(particle.getX())/pybeam.units.kpc,"kpc")
print("Time Z ->", particle.Time.z())

del particle # to be sure for flush tracer

fl = np.loadtxt("/tmp/trajectory.dat")
Rs = (fl[:,0]**2+fl[:,1]**2+fl[:,2]**2)**0.5
print("Rs.min() = ",Rs.min())

print("Max treads:",pybeam.omp_get_max_threads())

import os
import psutil

def getMemUsage():
	return psutil.Process(os.getpid()).memory_info().rss

#oldMem = getMemUsage()
#print("Mem usage:",oldMem)
#mf3 = pybeam.SavedMF("/tmp/magneticField.out")
#newMem = getMemUsage()
#print("Mem usage after one MF upload:",newMem," diff:",newMem-oldMem)
#for i in range(100):
#	mf3 = pybeam.SavedMF("/tmp/magneticField.out")
#
#newMem = getMemUsage()
#print("Mem usage after 100 MF upload:",newMem," diff:",newMem-oldMem)

def getValueGauss(x,y,z,t):
	B = np.array([x,y,z])+t.t()
	return list(B)

def getScale(t):
	return 0.1*pybeam.units.Mpc

custom_mf = pybeam.CustomMF(getValueGauss,getScale)

print("Call custom MF:", custom_mf.GetValueGauss([0,0,1],0), custom_mf.MinVariabilityScale(pybeam.CosmoTime(0))/pybeam.units.Mpc)

print("")
print("TESTS SUCCESSFULLY PASSED!")
