#!/bin/bash
# 16 random depth-3 sko49 nodes through dnntest, 3 at a time (the root run holds the 4th core),
# PRSM (-g 1.0 -G 1.0) + degree-4 filter (-d 4), 1 thread, run to 1,500 iterations or convergence
# (-H 0: no give-up), best known 23386 as the incumbent (cutoff 23384.05).
cd "$(dirname "$0")/../../../latest" || exit 1
OUT=../experiments/sko49_feasibility/depth3
P=${1:-3}
cat -n $OUT/nodes.txt | xargs -P $P -L 1 bash -c 'k=$0; spec=$1; stdbuf -oL ./dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 23386 -v -x $spec ../qaplib/sko49.dat > '$OUT'/node$k.log 2>&1'
echo "=== all done ==="
