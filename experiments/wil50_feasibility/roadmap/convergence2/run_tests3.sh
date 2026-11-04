#!/bin/bash
# convergence2, core 1 after run_tests.sh's core-1 E children (run_tests.sh's
# C n13, C n9, P 10 and P 2 on core 1 are skipped through placeholder logs,
# which this script removes before it runs C n9 and C n13 itself):
#  E2 child 49:45 of node 10, solved with node 10's equip pool, then the
#     search's cut phase from that start (inherited cuts at their multipliers,
#     300-iteration rounds, stall 3, quiet round of 1,500 from the best round;
#     the search runs the quiet round only when that round's primal is above
#     the cutoff, dnncuts always: read the round lines). In the search
#     (coord_search/trace.txt) this child ended 186 short, outside the zone,
#     and its subtree cost 137 nodes and 8,002 CPU-s.
#  then test C (Y12q) on n9 and n13, as in run_tests.sh.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
RUNS=$D/runs; CUTS=$R/experiments/dnn_cuts/dnncuts; W=$R/qaplib/wil50.dat
until grep -q "done E 47" "$RUNS/queue.txt" 2>/dev/null; do sleep 15; done
# run_tests.sh's core-1 queue passes the placeholders at once
sleep 15
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
f=$RUNS/E2_c49_phase.log
env DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 taskset -c 1 \
  "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -vc -I 4000 -H 300 -r 12 -i 300 -stall 3 -nocap -delta 50 \
  -freeze 0:1500 -loadpool "$RUNS/E_n10.pool" -x 39:27,30:48,25:49,49:45 "$W" > "$f" 2>&1
echo "$(date -u +%FT%TZ) done E2 49" >> "$RUNS/queue.txt"
rm -f "$RUNS/C_n9_Y12q.log" "$RUNS/C_n13_Y12q.log" "$RUNS/P_10_J1000.log" "$RUNS/P_2_J1000.log"
for spec in "n9 26:27,8:38,34:28" "n13 21:28,8:41,28:19"; do
  set -- $spec
  taskset -c 1 "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -delta 50 -nocap \
    -x $2 -r 12 -i 100 -freeze 0:500 -warmY "$W" > "$RUNS/C_$1_Y12q.log" 2>&1
  echo "$(date -u +%FT%TZ) done C $1 (core 1, run_tests3)" >> "$RUNS/queue.txt"
done
