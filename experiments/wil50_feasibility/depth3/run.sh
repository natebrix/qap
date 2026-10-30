#!/bin/bash
# 16 random depth-3 wil50 nodes through dnntest, 8 at a time, best-known 48816 as the incumbent
cd /Users/nathan_brixius/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1
OUT=/tmp/claude-501/wil50_d3
python3 - > $OUT/nodes.txt <<'PY'
import random
random.seed(20261029)
for k in range(16):
    f = random.sample(range(50), 3); l = random.sample(range(50), 3)
    print(",".join(f"{a}:{b}" for a, b in zip(f, l)))
PY
cat -n $OUT/nodes.txt | while read k spec; do echo "$k $spec"; done | \
  xargs -P 8 -L 1 bash -c 'k=$0; spec=$1; /usr/bin/time -p ./dnntest -f -T 1 -H 300 -i 800 -u 48816 -x $spec ../qaplib/wil50.dat > /tmp/claude-501/wil50_d3/node$k.log 2>&1'
echo "=== all 16 done ==="
grep -h "fathoming tolerance" $OUT/node1.log
for k in $(seq 1 16); do printf "node %2d %-14s " $k "$(sed -n ${k}p $OUT/nodes.txt)"; grep -o "lb=[0-9.]* status=[a-z_]* iters=[0-9]* cpu=[0-9.]* ms_per_iter=[0-9.]*.*fathomed=[01]" $OUT/node$k.log | sed 's/primal=[^ ]* //; s/certs=.*fathomed/fathomed/'; done
