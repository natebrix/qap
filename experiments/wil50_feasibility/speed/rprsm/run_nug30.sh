#!/bin/bash
# nug30 depth-2 nodes (m = 28) with the nug30 record setting's DNN switches
# (RUNNING.md: PRSM 1.0/1.0, degree 3, -L 2, beta 0.21, refresh 120, one
# CholQR pass, give-up at 200; cutoff from ub 6124), cold start, one node at a
# time on core 2.
#   usage: run_nug30.sh BINARY NAME:RPRSM [NAME:RPRSM ...]
# Logs: runs30/<name>_<id>.log, id = position in NODES30 below.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=$1; shift
# 4:10,6:8 and the two misses are the eigensolver README's PRSM nodes; the
# others are from experiments/dnn_bound/predictor/nug30_depth2.csv
NODES30=${NODES30:-"4:10,6:8 4:29,6:7 4:7,6:14 0:2,23:23 0:14,6:10 0:13,6:10 0:2,23:11 0:1,23:4"}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
mkdir -p "$D/runs30"
for cfg in "$@"; do
  name=${cfg%%:*}; rp=${cfg#*:}; [ "$rp" = "$cfg" ] && rp=0
  k=0
  for spec in $NODES30; do
    k=$((k + 1))
    f=$D/runs30/${name}_$k.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    if [ "$rp" != 0 ]; then export DNN_RPRSM=$rp; else unset DNN_RPRSM; fi
    taskset -c 2 "$BIN" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -B 0.21 -r 120 -H 200 -i 2500 -u 6124 -v \
      -x "$spec" "$R/qaplib/nug30.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
    echo "$name $k $spec $(grep RESULT "$f" | sed 's/.*lb=\([0-9.]*\).*status=\([a-z_]*\) iters=\([0-9]*\) cpu=\([0-9.]*\) ms_per_iter=\([0-9.]*\).*/\1 \2 \3 \4 \5/')"
  done
done
