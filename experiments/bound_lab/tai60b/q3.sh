#!/bin/bash
# q3.sh: second slot after q2: depth 4 under the on-path node (children 1:x of 44:45,43:58,45:26). First solves of
# 5 random children and 9 of the 13 same-cluster ones (facilities of the 14-point cluster in the best known
# solution; nearopt.py's excesses 0.25M-2.9M), then phases of four of those
D=$(cd $(dirname $0); pwd); cd $D
while pgrep -f "q2.sh" > /dev/null; do sleep 20; done
S=runs/d3_45_26path/out.state
./pass1.sh d4p1_under_path $S 4 5 4 1:39 1:59 1:24 1:29 1:36 1:35 1:32 1:0 1:23
for v in 39 59 29 35; do ./kid.sh d4_1_$v kidW $S 4 1:$v; done
