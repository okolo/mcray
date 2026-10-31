#!/bin/bash
if [ $# -gt 0 ]
then
    echo "usage: ./convert.sh <lowerdataIII.txt >lowerdataIII.dat" >&2
    exit 1
fi
#printing Z-scale in the format
# z1 z2 .. zMax
COUNTER=0
while [  $COUNTER -lt 401 ]; do
    echo $COUNTER | awk ' { printf("\t%.2f",0.02*$1) }'
    let COUNTER=COUNTER+1
done

echo

#printing data in the format
# E1	 dN/dE(z1,E1)	dN/dE(z2,E1)	dN/dE(z3,E1)	...			dN/dE(zMax,E1)
# E2	 dN/dE(z1,E2)	dN/dE(z2,E2)	dN/dE(z3,E2)	...			dN/dE(zMax,E2)
# ................
# Emax dN/dE(z1,Emax)	dN/dE(z2,Emax)	dN/dE(z3,Emax)	...			dN/dE(zMax,Emax)

awk ' { En=10^(0.19+0.01*NR); printf("%g\t", En); for (i=1; i<=NF; i++) printf("%g\t",1e-6/En*$i); printf("\n") }'
