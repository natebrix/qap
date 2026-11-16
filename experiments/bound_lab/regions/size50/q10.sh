#!/bin/bash
# after 34:36: do the children of a depth-4 survivor die at once? first solves of 10 random depth-5 children of
# 34:28 (branched 546k short); then the on-path depth-4 node 34:6 run long (3,000 phase iterations, no stall rule,
# wc = 2.88, incumbent best known + 1)
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "kid.sh d4_34_36" > /dev/null; do sleep 3; done
./pass1.sh d5p1_under_34_28 runs/d4_34_28/out.state 5 10 51
./kid.sh d4_34_6path_long kidLong runs/d3_32_0path/out.state 4 34:6
echo Q10 DONE
