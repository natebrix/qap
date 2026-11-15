#!/bin/bash
# queue 3: warm children on tai30b. Parent = depth-2 node, cold, 600 iterations; three depth-3 children
# warm (Z and conditioned Y, -Q b), 600 iterations, no give-up; the warm-call beta factor w and dnn.dscale vary.
# P1 = 0:3,1:7 is on the optimal path (QAPSOLN/tai30b.sln: facility i at location p[i]); child 2:10 continues it.
cd "$(dirname "$0")"
WT=/home/user/wt-admm/latest/dnntest
for P in "P1 0:3,1:7 2:10;2:14;5:6" "P2 4:9,11:20 2:10;7:3;15:15"; do
  for cfg in w0.36 w1 w1_ds0.5 w2.4 w1.67_ds0.5 w0.36_ds0.5; do
    set -- $P
    [ -s logs/tai30b_$1_$cfg.log ] && continue
    $WT -p par_$cfg.par -T 1 -x $2 -J 600 -Q b -H 0 -i 600 -k 10 -u 637117113 -v -K "$3" ../../../qaplib/tai30b.dat > logs/tai30b_$1_$cfg.log 2>&1
  done
done
