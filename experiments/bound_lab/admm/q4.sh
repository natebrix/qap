#!/bin/bash
# queue 4: tai50b, warm child of depth-3 node 7 (NOTES.md section 2: 0:44,49:49,25:31). Parent cold 600 iterations,
# child (depth 4, m = 46) warm (Z and conditioned Y), 600 iterations, no give-up: production settings against dnn.dscale = 0.5
cd "$(dirname "$0")"
while pgrep -f "extrap.py|q3.sh" > /dev/null; do sleep 5; done
WT=/home/user/wt-admm/latest/dnntest
for cfg in w0.36 w0.36_ds0.5; do
  $WT -p par_$cfg.par -T 1 -x 0:44,49:49,25:31 -J 600 -Q b -H 0 -i 600 -k 10 -u 458821517 -v -K "1:23" ../../../qaplib/tai50b.dat > logs/tai50b_n7kid_$cfg.log 2>&1
done
./q2.sh
