#!/usr/bin/env python3
"""Independently integrate prepared Fe56 cross sections on an analytic CMB.

This diagnostic uses composite Simpson quadrature, not mcray's background
integrals or GSL. It reports differences; it does not define a TALYS physical
acceptance criterion. Constants deliberately match the existing mcray Units.
"""
import argparse
import math
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--scan',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
p=args.data
header=args.scan.read_text().splitlines()[0]
if 'T_K=2.7255 ' not in header or 'background=cmb ' not in header or 'primary=A56Z26 ' not in header:
 raise ValueError('Expected a frozen 2.7255 K CMB scan of Fe56')
grid=[float(l) for l in (p/'grid_talys.tsv').read_text().splitlines() if not l.startswith('#')]
row=next(l.split() for l in (p/'total.tsv').read_text().splitlines() if l.startswith('56 26 '));sigma=list(map(float,row[2:]));mass=float(next(l.split()[2] for l in (p/'masses.tsv').read_text().splitlines() if l.startswith('56 26 ')))
T=2.7255/1.16e10;hbarc=6.582122020e-22*2.99792458e10;mpc=3.0856775807e24

def independent(gamma,n):
 integral=0
 for k0,k1,s0,s1 in zip(grid,grid[1:],sigma,sigma[1:]):
  if s0==s1==0:continue
  dl=math.log(k1/k0);total=0
  for j in range(n+1):
   w=j/n;k=k0*math.exp(dl*w)
   s=math.exp(math.log(s0)*(1-w)+math.log(s1)*w) if s0>0 and s1>0 else s0*(1-w)+s1*w
   x=k/(2*gamma*T)
   photon=0 if x>700 else (-math.log1p(-math.exp(-x)) if x>1 else -math.log(-math.expm1(-x)))
   total+=(1 if j in (0,n) else 4 if j%2 else 2)*k*k*s*photon
  integral+=total*dl/(3*n)
 return integral *1e-27*T/(2*math.pi**2*gamma*gamma*hbarc**3)*mpc
results=[]
for line in args.scan.read_text().splitlines():
 if line.startswith('#'):continue
 energy,rate,_=map(float,line.split()[:3]);gamma=energy*56e-6/mass
 if rate<1e-10:continue
 a,b=independent(gamma,128),independent(gamma,256)
 results.append((energy,rate,b,rate/b-1,abs(a/b-1)))
if not results: raise ValueError('No measurable rates in scan')
print('points',len(results),'max relative difference',max(abs(r[3]) for r in results),'reference quadrature convergence',max(r[4] for r in results))
print(*results[::10],sep='\n')
args.output.write_text('# E_per_nucleon_eV mcray_rate_Mpc^-1 independent_rate_Mpc^-1 rel_difference quadrature_change\n'+'\n'.join(' '.join(format(x,'.17g') for x in r) for r in results)+'\n')
