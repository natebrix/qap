#!/bin/bash
# ./mkpar.sh NAME INSTANCE REGIONFAMILY [extra "key = value" lines]: par/NAME.par = tai50b_size/par/A_wc.par
# (tai50b.par + the cut block + dnn.beta_sched = w=0.36,wc=2.88) with the instance's depth range and region file
H=$(cd "$(dirname "$0")" && pwd); N=$1; I=$2; F=$3; shift 3; n=${I//[a-z]/}
sed -e "s/^dnn.depths = .*/dnn.depths = 0:$((n-2))/" -e "s#^cut.regions = .*#cut.regions = $H/regions/$I.$F.txt#" $H/../tai50b_size/par/A_wc.par > $H/par/$N.par
for l in "$@"; do k=${l%% =*}; grep -q "^$k = " $H/par/$N.par && sed -i "s#^$k = .*#$l#" $H/par/$N.par || echo "$l" >> $H/par/$N.par; done
