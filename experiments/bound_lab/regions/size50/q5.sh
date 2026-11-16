#!/bin/bash
# after q4: the furthest random depth-1 child (31:16, 7.49M short) under wc = 2.88, cap 1,000
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "q4.sh" > /dev/null; do sleep 5; done
./kid.sh d1_31_16_wc288 kidW288 runs/rootL/root.state 1 31:16
echo Q5 DONE
