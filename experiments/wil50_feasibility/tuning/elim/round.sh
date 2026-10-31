#!/bin/bash
# $1 = node id, $2 = spec. Round 1: DNN bound with child bounds; round 2: re-bound with the children the certificate fathoms eliminated.
S=${S:-$(dirname "$0")}
cd /Users/nathan_brixius/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1
k=$1; spec=$2
./dnntest -T 1 -d 4 -g 1.0 -G 1.0 -H 0 -i 2500 -u 48816 -C 2 -x $spec ../qaplib/wil50.dat > $S/r1_$k.log 2>&1
el=$(awk -v c=48814.05 '$1=="CHILD" && $4>c {printf "%s%d:%d", (n++?",":""), $2, $3}' $S/r1_$k.log)
echo "$el" > $S/elim_$k.txt
if [ -n "$el" ]; then
  ./dnntest -T 1 -d 4 -g 1.0 -G 1.0 -H 0 -i 2500 -u 48816 -C 2 -E "$el" -x $spec ../qaplib/wil50.dat > $S/r2_$k.log 2>&1
  el2=$(awk -v c=48814.05 '$1=="CHILD" && $4>c {n++} END{print n+0}' $S/r2_$k.log)
else
  el2=0
fi
r1=$(grep -o 'lb=[0-9.]* primal=[0-9.]* status=[a-z_]* iters=[0-9]* cpu=[0-9.]*' $S/r1_$k.log)
c1=$(grep -o 'CHILDREN.*' $S/r1_$k.log)
r2=$(grep -o 'lb=[0-9.]* primal=[0-9.]* status=[a-z_]* iters=[0-9]* cpu=[0-9.]*' $S/r2_$k.log 2>/dev/null)
echo "node $k $spec | R1 $r1 | $c1 | R2(elim $(echo $el | tr ',' '\n' | grep -c .)) $r2 | R2 children over cutoff: $el2"
