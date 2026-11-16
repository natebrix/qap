#!/bin/bash
# depth-4 phases (wc = 2.88, cap 1,000) under the on-path depth-3 node: the furthest, the nearest and a middle
# one of the 8 random children that did not fathom at the first solve
cd /home/user/qap/experiments/bound_lab/regions/size50
./kid.sh d4_34_28 kidW288 runs/d3_32_0path/out.state 4 34:28
./kid.sh d4_34_19 kidW288 runs/d3_32_0path/out.state 4 34:19
./kid.sh d4_34_36 kidW288 runs/d3_32_0path/out.state 4 34:36
echo Q8 DONE
