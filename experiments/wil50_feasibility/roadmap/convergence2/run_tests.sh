#!/bin/bash
# convergence2 tests (REPORT.md section 4). Two single-threaded queues pinned to
# cores 0 and 1 (the compute rules of the brief: cores 2-3 belong to another
# agent). Waits until the node-10 coordinator measurement has finished.
#
#  C  cut phase with short warm-Y rounds: the budget study's protocol
#     (same_node_continuation/budget_study: -I 400, 1,000 cuts a round, delta 50
#     uncapped, dnncuts' default eigensolver settings) with -warmY, 12 rounds of
#     100 iterations and a quiet round of 500 ("Y12q"); compared with that
#     study's B / Y12 / Y12s logs (round 0 must match theirs to the decimal).
#  E  a cheap cut phase at the clear miss node 10 (the search's protocol: 300
#     uncut iterations, 300-iteration rounds, stall 3, no quiet round) that
#     only equips its children with a pool; then 13 of its depth-4 children
#     (the 9 that missed in coord_search/trace.txt and 4 that fathomed) solve
#     once with that pool frozen, cold, give-up after 300 as in the search.
#  S  a warm start shipped at reduced precision: parent 13 runs exactly 300
#     iterations (the state warm_y/runs used: its give-up fired at 300), its
#     saved Y and Z are rounded to bf16, then its 4 children start from the
#     conditioned Y; compared with warm_y/runs child_5..8 (start b).
#  P  parent investment: depth-3 parents 13, 7, 10, 2 run 1,000 iterations
#     (no give-up) instead of giving up at 300, then their 4 depth-4 children
#     each (stronger_bounds/depth4/nodes.txt) from the conditioned parent Y;
#     compared with warm_y/runs (the same children after a parent at 300).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
RUNS=$D/runs; mkdir -p "$RUNS"
CUTS=$R/experiments/dnn_cuts/dnncuts; DT=$R/latest/dnntest; W=$R/qaplib/wil50.dat
GATE=/home/user/qap/experiments/wil50_feasibility/coord_search/coordinator_stacked.log
until grep -q "coordinator: 4 slots" "$GATE" 2>/dev/null; do sleep 60; done
echo "$(date -u +%FT%TZ) start" >> "$RUNS/queue.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PROD="DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40"

C() {  # core, node id, fixes
  local f=$RUNS/C_$2_Y12q.log
  [ -s "$f" ] && grep -q FINAL "$f" && return
  taskset -c $1 "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -delta 50 -nocap \
    -x $3 -r 12 -i 100 -freeze 0:500 -warmY "$W" > "$f" 2>&1
  echo "$(date -u +%FT%TZ) done C $2" >> "$RUNS/queue.txt"
}
E0() {
  local f=$RUNS/E_n10_phase.log
  [ -s "$RUNS/E_n10.pool" ] && return
  env $PROD taskset -c $1 "$CUTS" -grid 5:10 -u 48816 -T 1 -k 1000 -v -vc -I 300 -delta 50 -nocap \
    -stall 3 -x 39:27,30:48,25:49 -r 12 -i 300 -savepool "$RUNS/E_n10.pool.tmp" "$W" > "$f" 2>&1
  mv "$RUNS/E_n10.pool.tmp" "$RUNS/E_n10.pool"
  echo "$(date -u +%FT%TZ) done E0" >> "$RUNS/queue.txt"
}
E() {  # core, child facility (placed on location 45)
  until [ -s "$RUNS/E_n10.pool" ]; do sleep 30; done
  local f=$RUNS/E_c$2.log
  [ -s "$f" ] && grep -q FINAL "$f" && return
  env $PROD taskset -c $1 "$CUTS" -grid 5:10 -u 48816 -T 1 -v -vc -I 4000 -H 300 -r 0 \
    -loadpool "$RUNS/E_n10.pool" -x 39:27,30:48,25:49,$2:45 "$W" > "$f" 2>&1
  echo "$(date -u +%FT%TZ) done E $2" >> "$RUNS/queue.txt"
}
P() {  # core, parent id, parent fixes, children
  local f=$RUNS/P_$2_J1000.log
  [ -s "$f" ] && [ "$(grep -c '^RESULT' "$f")" -ge 5 ] && return
  env $PROD taskset -c $1 "$DT" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 1500 -u 48816 -v \
    -x $3 -J 1000 -K $4 -Q b "$W" > "$f" 2>&1
  echo "$(date -u +%FT%TZ) done P $2" >> "$RUNS/queue.txt"
}

S() {  # core: parent 13 at exactly 300 iterations, state rounded to bf16, its 4 children
  local f=$RUNS/S_13_J300_bf16.log
  [ -s "$f" ] && [ "$(grep -c '^RESULT' "$f")" -ge 5 ] && return
  env DNN_AMX=1 DNN_CHOLQR_PASSES=1 taskset -c $1 "$DT" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 1500 \
    -u 48816 -v -x 21:28,8:41,28:19 -J 300 -Wr 16 -K 3:2,44:2,15:2,18:2 -Q b "$W" > "$f" 2>&1
  echo "$(date -u +%FT%TZ) done S 13" >> "$RUNS/queue.txt"
}

q0() {
  E0 0
  for c in 48 49 38 29 28 19 18; do E 0 $c; done
  C 0 n14 47:23,10:44,40:34
  C 0 n6 20:12,40:28,17:25
  S 0
  C 0 n10 39:27,30:48,25:49
}
q1() {
  P 1 13 21:28,8:41,28:19 3:2,44:2,15:2,18:2
  P 1 7 12:7,33:6,32:11 5:34,8:34,42:34,0:34
  for c in 9 8 0 46 37 47; do E 1 $c; done   # 9, 8 missed; 0, 46, 37, 47 fathomed at 130/200/290/450
  C 1 n13 21:28,8:41,28:19
  C 1 n9 26:27,8:38,34:28
  P 1 10 39:27,30:48,25:49 7:5,26:5,3:5,47:5
  P 1 2 18:43,11:29,46:44 8:21,9:21,30:21,37:21
}
q0 & q1 & wait
echo "$(date -u +%FT%TZ) all done" >> "$RUNS/queue.txt"
