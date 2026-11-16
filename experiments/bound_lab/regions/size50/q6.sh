#!/bin/bash
# after q5: the on-path depth-3 node (incumbent best known + 1), then the first solves of 10 random depth-4
# children (+ on-path) below it
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "q5.sh" > /dev/null; do sleep 5; done
./kid.sh d3_32_0path kidP runs/d2_33_32path/out.state 3 32:0
BR=$(awk '{print $10}' runs/d3_32_0path/nodes.txt)
echo "depth-3 on-path node branches on $BR"
L=${BR#?}; P=$(python3 -c "
t=open('/home/user/qap/QAPSOLN/tai50b.sln').read().replace(',',' ').split(); p=[int(x)-1 for x in t[2:52]]
print('$L:%d' % p[$L] if '$BR'[0]=='r' else '%d:$L' % p.index($L))")
./pass1.sh d4p1_under_path runs/d3_32_0path/out.state 4 10 41 $P
echo Q6 DONE
