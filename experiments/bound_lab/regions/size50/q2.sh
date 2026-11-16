#!/bin/bash
# depth-2 first solves: 10 random children (+ the on-path one) of the on-path depth-1 node, 10 of the median one
cd /home/user/qap/experiments/bound_lab/regions/size50
./pass1.sh d2p1_under_31_7 runs/d1_31_7path/out.state 2 10 21 33:32
./pass1.sh d2p1_under_31_23 runs/d1_31_23/out.state 2 10 22
echo Q2 DONE
