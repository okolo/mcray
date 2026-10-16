#!/usr/bin/env python3
"""Independent-history standard errors; bins are descriptive, not acceptance tests."""
import argparse
import math
import pathlib
import re
p=argparse.ArgumentParser()
p.add_argument('old',type=pathlib.Path)
p.add_argument('current',type=pathlib.Path)
a=p.parse_args()
def read(path):
    text=path.read_text()
    n=int(re.search(r' N=(\d+)',text)[1])
    rows=[list(map(float,line.split())) for line in text.splitlines() if not line.startswith('#')]
    if len(rows)!=32 or n<2:raise ValueError('Incomplete ensemble')
    metadata=dict(re.findall(r"(\w+)=([^\s]+)",text.splitlines()[0]))
    if any(not math.isfinite(v) for row in rows for v in row):raise ValueError("Non-finite output")
    return n,rows,metadata
def stats(n,s,ss):
    return s/n, math.sqrt(max(0.,(ss-s*s/n)/(n-1)/n))
n,old,old_meta=read(a.old)
m,new,new_meta=read(a.current)
for key in ("E0_eV","Emin_eV","L_Mpc","source_z","CMB_K","tolerance","max_rate_variation"):
    if key not in old_meta or key not in new_meta or old_meta[key]!=new_meta[key]:raise ValueError("Scenario mismatch: "+key)
if old_meta["seed"]==new_meta["seed"]:raise ValueError("Use independent seeds for independent-ensemble errors")
print('# species log10_E_low old_count current_count current_minus_old SE_difference z_score old_energy_eV current_energy_eV energy_SE_difference energy_z_score')
for l,r in zip(old,new):
    if l[0]!=r[0]:raise ValueError('Bin mismatch')
    index=int(l[0]);x,sx=stats(n,l[1],l[2]);y,sy=stats(m,r[1],r[2]);se=math.hypot(sx,sy)
    if x+y==0:continue
    ex,sex=stats(n,l[3],l[4]);ey,sey=stats(m,r[3],r[4]);ese=math.hypot(sex,sey)
    print('photon' if index<16 else 'lepton',12+(index%16)/4,x,y,y-x,se,(y-x)/se if se else 'NA',ex,ey,ese,(ey-ex)/ese if ese else 'NA')
