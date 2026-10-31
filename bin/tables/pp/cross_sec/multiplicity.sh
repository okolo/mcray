#!/bin/bash
cat nu_pp | awk -f integral_sigma_nu.awk >sigma_tot_nu
for E in `cat sigma_in | awk 'NF==2 { print $1 }'`
do
    E_in_eV=`echo $E | awk '{printf("%g",$1*1e9)}'`
    F=`cat sigma_in | python func.py $E logscale`
    nF=`cat sigma_tot_nu | python func.py $E_in_eV logscale`
    echo $E $F $nF | awk '{printf("%g\t%g\n",$1,$3/$2)}'
done
rm sigma_tot_nu
