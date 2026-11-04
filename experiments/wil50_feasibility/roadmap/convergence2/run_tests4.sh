#!/bin/bash
# convergence2, core 1 (replaces run_tests3.sh, stopped at 22:57 during E2's quiet
# round, which the search would not run: the best round's primal 48,767.8 is
# below the cutoff):
#  E2w child 49:45 of node 10 with node 10's equip pool, then the search's phase
#      from that start as in E2 but with -warmY (each round and the quiet round
#      restart from the previous Y as well as Z, QAP_DNN_CUT_WARMY=1 in bnb.C).
#  then test C (Y12q) on n9 and n13, as in run_tests.sh.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
RUNS=$D/runs; CUTS=$R/experiments/dnn_cuts/dnncuts; W=$R/qaplib/wil50.dat
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
f=$RUNS/E2w_c49_phase_warmY.log
env DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 taskset -c 1 \
  "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -vc -I 4000 -H 300 -r 12 -i 300 -stall 3 -nocap -delta 50 \
  -freeze 0:1500 -warmY -loadpool "$RUNS/E_n10.pool" -x 39:27,30:48,25:49,49:45 "$W" > "$f" 2>&1
echo "$(date -u +%FT%TZ) done E2w 49" >> "$RUNS/queue.txt"
rm -f "$RUNS/C_n9_Y12q.log" "$RUNS/C_n13_Y12q.log" "$RUNS/P_10_J1000.log" "$RUNS/P_2_J1000.log"
for spec in "n9 26:27,8:38,34:28" "n13 21:28,8:41,28:19"; do
  set -- $spec
  taskset -c 1 "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -delta 50 -nocap \
    -x $2 -r 12 -i 100 -freeze 0:500 -warmY "$W" > "$RUNS/C_$1_Y12q.log" 2>&1
  echo "$(date -u +%FT%TZ) done C $1 (core 1)" >> "$RUNS/queue.txt"
done
