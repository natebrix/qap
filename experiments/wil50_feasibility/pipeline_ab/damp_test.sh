#!/bin/bash
# First-step damping (QAP_DNN_CUT_DAMP1, dnncuts -damp1) on the case that
# motivated it: child 49:45 of node 10 with node 10's equip-phase pool
# (roadmap/convergence2 E2/E2w: the first Polyak step from the inherited
# multipliers drops the bound 88-116 units). The warmy arm's phase settings:
# Y kept between rounds, 12 rounds of 200 iterations, quiet round at most
# 1,000, stall rule 3, uncapped target. Two runs, damp1 1.0 and 0.3, one core
# each. Logs in damp/.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
CUTS=$R/experiments/dnn_cuts/dnncuts; W=$R/qaplib/wil50.dat
POOL=$R/experiments/wil50_feasibility/roadmap/convergence2/runs/E_n10.pool
mkdir -p "$D/damp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export DNN_AMX=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
for f in ${DAMPS:-1.0 0.3}; do
  "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 4000 -H 300 -r 12 -i 200 -stall 3 -nocap -delta 50 \
    -freeze 0:1000 -warmY -damp1 "$f" -loadpool "$POOL" -x 39:27,30:48,25:49,49:45 "$W" \
    > "$D/damp/c49_damp${f/:/_}.log" 2>&1 &
done
wait
grep -h -E "^  round|FINAL|stalled|^quiet|fathom" "$D"/damp/c49_damp*.log
