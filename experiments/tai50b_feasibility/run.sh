#!/bin/bash
# tai50b feasibility runs on the 2-core cloud Xeon (OpenBLAS, 1 thread per process).
# Root and 12 random depth-3 nodes through dnntest with the sko49 note's settings:
# PRSM (-g 1.0 -G 1.0), degree-4 filter (-d 4), 1 thread, no give-up (-H 0),
# 1,500 iterations or convergence, best known 458,821,517 as the incumbent.
# The root takes one core; the depth-3 nodes queue on the other, then on both.
cd "$(dirname "$0")/../../latest" || exit 1
OUT=../experiments/tai50b_feasibility
UB=458821517
export OPENBLAS_NUM_THREADS=1
DNN="./dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u $UB -v"
run_node() { k=$1; spec=$(sed -n "${k}p" $OUT/depth3/nodes.txt); stdbuf -oL $DNN -x $spec ../qaplib/tai50b.dat > $OUT/depth3/node$k.log 2>&1; }
# core 1: the root, then nodes 7-12; core 2: nodes 1-6
( stdbuf -oL $DNN ../qaplib/tai50b.dat > $OUT/root/root.log 2>&1; for k in 7 8 9 10 11 12; do run_node $k; done ) &
( for k in 1 2 3 4 5 6; do run_node $k; done ) &
wait
echo "=== all done ==="
