#!/bin/bash
# depth-1 phases (cap 1,000): the nearest child, a median one, the on-path one (incumbent best known + 1)
cd /home/user/qap/experiments/bound_lab/regions/size50
./kid.sh d1_31_1 kidP runs/rootL/root.state 1 31:1
./kid.sh d1_31_7path kidP runs/rootL/root.state 1 31:7 458821518
./kid.sh d1_31_23 kidP runs/rootL/root.state 1 31:23
echo Q1 DONE
