#!/bin/bash
dlogE=`cat nu_pp | awk 'BEGIN{prev=1}NR<3{print $2/prev; prev=$2}' | awk 'NR==2 {print log($1)}'`
echo dlogE=$dlogE
Elast=`tac nu_pp | awk 'NR==5 {print $1}'`
echo Elast=$Elast
sigmaDiffLast=`cat nu_pp | awk 'BEGIN{sum=0} $1==Elast{ sum=sum+$3+$4+$5+$6 }END{print sum*dlogE}' Elast=$Elast dlogE=$dlogE`
echo sigmaDiffLast=$sigmaDiffLast
ElastGeV=`echo $Elast | awk '{print 1e-9*$1}'`
echo ElastGeV=$ElastGeV
sigmaTotLast=`cat sigma_in | python func.py $ElastGeV logscale`
echo sigmaTotLast=$sigmaTotLast
echo $sigmaTotLast $sigmaDiffLast | awk '{print $2/$1}'
