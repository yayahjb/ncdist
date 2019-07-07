#!/bin/bash
latsym1=$1
a1=$2
b1=$3
c1=$4
alpha1=$5
beta1=$6
gamma1=$7
shift 7
latsym2=$1
a2=$2
b2=$3
c2=$4
alpha2=$5
beta2=$6
gamma2=$7
s6dist_val=`./s6dist_app $latsym1 $a1 $b1 $c1 $alpha1 $beta1 $gamma1 $latsym2 $a2 $b2 $c2 $alpha2 $beta2 $gamma2 | tail -1`
cs6dist_val=`./cs6dist_app $latsym1 $a1 $b1 $c1 $alpha1 $beta1 $gamma1 $latsym2 $a2 $b2 $c2 $alpha2 $beta2 $gamma2 | tail -1`
ncdist_val=`./ncdist $latsym1 $a1 $b1 $c1 $alpha1 $beta1 $gamma1 $latsym2 $a2 $b2 $c2 $alpha2 $beta2 $gamma2 | tail -1`
#echo $s6dist_val
#echo $cs6dist_val
#echo $ncdist_val
echo "cell1: [ $latsym1 $a1 $b1 $c1 $alpha1 $beta1 $gamma1] cell2: [ $latsym2 $a2 $b2 $c2 $alpha2 $beta2 $gamma2] s6dist: $s6dist_val cs6dist: $cs6dist_val ncdist: $ncdist_val"
