#!/bin/bash
# the phase penalty (dnn.beta_sched wc, cpfast/README.md) on the median depth-1 child 31:23, warm from the root:
# baseline (no wc, the phase at w = 0.36): 359,412 short after 1,000 phase iterations
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f q2.sh > /dev/null; do sleep 5; done
./kid.sh d1_31_23_wc288 kidW288 runs/rootL/root.state 1 31:23
./kid.sh d1_31_23_wc100 kidW10 runs/rootL/root.state 1 31:23
echo Q3 DONE
