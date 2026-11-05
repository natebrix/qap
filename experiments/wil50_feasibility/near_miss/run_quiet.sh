#!/bin/bash
# The near miss [0:25 9:21 49:45] (depth 6; the cloud stacked run's phase
# ended 0.27 short after its 1,500-iteration quiet round, the laptop base's
# 9.6 short): the production phase protocol from a cold first solve, with a
# quiet round of 3,000 iterations and the certified bound at every check
# (-vc), to see whether and when a longer quiet round crosses the cutoff.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
taskset -c ${CORE:-3} "$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -u 48816 -T 1 -k 1000 -v -vc -I 4000 -H 300 \
  -r 12 -i 300 -stall 3 -nocap -delta 50 -freeze 0:3000 -x 39:27,30:48,25:49,49:45,0:25,9:21 \
  "$R/qaplib/wil50.dat" > "$D/logs/n0-25_9-21_quiet3000.log" 2>&1
