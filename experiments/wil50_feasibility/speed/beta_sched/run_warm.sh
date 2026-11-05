#!/bin/bash
# Warm-started children: the wil50 parent 26:27,8:38,34:28 (all32 node 9, a
# miss) runs as in run_all32.sh, then four of its children (all32 nodes 17,
# 18, 19, 20: 37:14, 37:37, 37:21, 37:26) start from its final state as in
# the search with QAP_DNN_WARM_Y=2 (dnntest -Q b: Z and the conditioned Y).
# Same flags as run_all32.sh, one process on core 1 per configuration.
# CONFIGS="name@env ..." as in run_all32.sh; DNN_BETA_SCHED's nodes=warm
# keeps the parent identical to the baseline's. PARENT, KIDS and SUFFIX select
# another family: PARENT=20:12,40:28,17:25 (all32 node 6) KIDS=43:18,43:7,43:22,43:19
# (all32 nodes 25-28) SUFFIX=_p6. Logs: warm/<name><SUFFIX>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CONFIGS=${CONFIGS:-"base@-"}
KIDS=${KIDS:-"37:14,37:37,37:21,37:26"}
PARENT=${PARENT:-"26:27,8:38,34:28"}
mkdir -p "$D/warm"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 ${envs//+/ } \
    taskset -c 1 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -u 48816 -v \
    -x "$PARENT" -K "$KIDS" -Q b "$R/qaplib/wil50.dat" > "$D/warm/$name$SUFFIX.log" 2>&1
done
