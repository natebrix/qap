#!/bin/bash
# depth-2 phases under the on-path depth-1 node 31:7: the two furthest random children and the on-path child;
# then the first solves of 10 random depth-3 children (+ on-path) of the on-path depth-2 node
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f q3.sh > /dev/null; do sleep 5; done
PARK=$(cat bestpar.txt 2>/dev/null || echo kidP)
./kid.sh d2_33_28 $PARK runs/d1_31_7path/out.state 2 33:28
./kid.sh d2_33_9 $PARK runs/d1_31_7path/out.state 2 33:9
./kid.sh d2_33_32path kidP runs/d1_31_7path/out.state 2 33:32
BR=$(awk '{print $10}' runs/d2_33_32path/nodes.txt)   # the on-path depth-2 node's branching line, e.g. r32
echo "depth-2 on-path node branches on $BR"
L=${BR#?}; P=$(python3 -c "
t=open('/home/user/qap/QAPSOLN/tai50b.sln').read().replace(',',' ').split(); p=[int(x)-1 for x in t[2:52]]
print('$L:%d' % p[$L] if '$BR'[0]=='r' else '%d:$L' % p.index($L))")
./pass1.sh d3p1_under_path runs/d2_33_32path/out.state 3 10 31 $P
echo Q4 DONE
