#!/bin/bash
# E1: certificate variants at equal iterates on fathoming first solves.
# Children of node 10 (39:27,30:48,25:49 + X:45), cold, as in
# pool_search/trace.txt and coord_search/trace.txt (warm 0 there), with the
# production solver flags (PRSM 1/1, degree 3, tight lambda 2, refresh 40, AMX
# bundle D, give-up 300). A child given as a list with ':' ("49:45,0:12") is
# taken as extra fixes under node 10 (a cold depth-5 node).
# DNN_CERTX evaluates, at the same iterates, the exact-LAP GL, Hahn-Grant
# skew polishing (8 passes), the reduce-pass polishing (12 + 1), the average
# of the last 4 multipliers, the projected post-step multiplier and (bit 4)
# the line search toward the post-step multiplier; the iterates are not
# changed (latest/dnn.C, "Experimental certificate variants").
# Usage: run_e1.sh <core> <mode> <child> [<child> ...]   (child = X for X:45)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; mode=$2; shift 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
for x in "$@"; do
  case $x in *:*) fx="39:27,30:48,25:49,$x";; *) fx="39:27,30:48,25:49,$x:45";; esac
  f=$D/logs/e1_c$(echo "$x" | tr ':,' '-_')_m${mode}.log
  grep -q '^RESULT' "$f" 2>/dev/null && continue
  echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core fixes $fx" > "$f.tmp"
  DNN_CERTX=$mode DNN_CERTX_EVERY=${EVERY:-5} DNN_CERTX_WINDOW=${WIN:-50} \
    taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 \
    -x "$fx" "$R/qaplib/wil50.dat" >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
done
