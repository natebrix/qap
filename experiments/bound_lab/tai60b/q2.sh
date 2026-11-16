#!/bin/bash
# q2.sh: second tai60b slot, after the long root: depth-2 first solves under the on-path node, then the phases
# of the two random depth-1 survivors of d1p1
D=$(cd $(dirname $0); pwd); cd $D
until grep -q "^exit" runs/rootL/cmd.txt; do sleep 20; done
./pass1.sh d2p1_under_path runs/d1_44_45path/out.state 2 12 2
./kid.sh d1_44_58 kidW runs/rootM/root.state 1 44:58
./kid.sh d1_44_39 kidW runs/rootM/root.state 1 44:39
