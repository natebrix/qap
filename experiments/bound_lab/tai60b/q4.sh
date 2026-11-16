#!/bin/bash
# q4.sh: second slot after q3: how fast do a survivor's children die? First solves of 10 random children of the
# depth-4 survivor 1:39 (same-cluster facility; branched 0.40M short), plus its own "on-path-like" child
D=$(cd $(dirname $0); pwd); cd $D
while pgrep -f "q3.sh" > /dev/null; do sleep 20; done
./pass1.sh d5p1_under_1_39 runs/d4_1_39/out.state 5 10 5
