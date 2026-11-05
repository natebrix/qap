#!/bin/bash
# A/B of DNN_RPRSM on the 32 wil50 depth-3 nodes of tuning/giveup/all32.txt,
# adapted from solvers/accel/run_all32.sh: production-like flags (degree 3,
# -L 2, PRSM -g 1.0 -G 1.0, -r 40, -H 300, -i 2500 -u 48816,
# DNN_CHOLQR_PASSES=1), one node at a time (-P 1) on core 2.
#   usage: run_all32.sh BINARY NAME:RPRSM [NAME:RPRSM ...]
# NODES="19 21 ..." restricts to those line numbers of all32.txt (default all).
# Configurations run one after the other, node by node. Resumable: a node whose
# log has a RESULT line is skipped. Logs: runs/<name>_<id>.log.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=$1; shift
LIST="$R/experiments/wil50_feasibility/tuning/giveup/all32.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
mkdir -p "$D/runs"
NODES=${NODES:-$(seq 1 32)}
for cfg in "$@"; do
  name=${cfg%%:*}; rp=${cfg#*:}; [ "$rp" = "$cfg" ] && rp=0
  for k in $NODES; do
    spec=$(sed -n "${k}p" "$LIST")
    f=$D/runs/${name}_$k.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    if [ "$rp" != 0 ]; then export DNN_RPRSM=$rp; else unset DNN_RPRSM; fi
    taskset -c 2 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -u 48816 -v \
      -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
    echo "$name $k $(grep RESULT "$f" | sed 's/.*status=\([a-z_]*\) iters=\([0-9]*\) cpu=\([0-9.]*\).*/\1 \2 \3/')"
  done
done
