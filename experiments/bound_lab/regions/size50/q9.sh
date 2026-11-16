#!/bin/bash
# a depth-4 child whose facility belongs to another cluster in the best known solution (34:36, 2.36M short), then
# the on-path depth-4 node 34:6 run long (3,000 phase iterations, no stall rule, wc = 2.88, incumbent best known + 1):
# where does the bound of a node on the path flatten?
cd /home/user/qap/experiments/bound_lab/regions/size50
./kid.sh d4_34_36 kidW288 runs/d3_32_0path/out.state 4 34:36
./kid.sh d4_34_6path_long kidLong runs/d3_32_0path/out.state 4 34:6
echo Q9 DONE
