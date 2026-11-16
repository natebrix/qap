#!/bin/bash
# after q6 (one solver process at a time): children of the far depth-1 child 31:16 (branched 1.39M short),
# first solves of 10 random; then the on-path depth-1 node again under wc = 2.88
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "q6.sh" > /dev/null; do sleep 5; done
./pass1.sh d2p1_under_31_16 runs/d1_31_16_wc288/out.state 2 10 23
./kid.sh d1_31_7path_wc288 kidW288 runs/rootL/root.state 1 31:7 458821518
echo Q7 DONE
