#!/bin/bash
# E1: certificate variants at equal iterates on fathoming first solves.
# Depth-4 children of node 10 (39:27,30:48,25:49 + X:45), cold, as in
# pool_search/trace.txt (warm 0 there), with the production solver flags
# (PRSM 1/1, degree 3, tight lambda 2, refresh 40, AMX bundle D, give-up 300).
# DNN_CERTX evaluates, at the same iterates, the exact-LAP GL, Hahn-Grant
# skew polishing (8 passes), the average of the last 4 multipliers and (mode
# 15) the line search toward the post-step multiplier; the iterates are not
# changed (latest/dnn.C, "Experimental certificate variants").
# Usage: run_e1.sh <core> <mode> <child> [<child> ...]   (child = X for X:45)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; mode=$2; shift 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
for x in "$@"; do
  f=$D/logs/e1_c${x}_m${mode}.log
  grep -q '^RESULT' "$f" 2>/dev/null && continue
  echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core" > "$f.tmp"
  DNN_CERTX=$mode DNN_CERTX_EVERY=5 DNN_CERTX_WINDOW=80 \
    taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 \
    -x "39:27,30:48,25:49,$x:45" "$R/qaplib/wil50.dat" >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
done
