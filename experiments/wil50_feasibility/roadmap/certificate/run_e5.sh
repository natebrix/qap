#!/bin/bash
# E5: does polishing longer pay? Child 37:45 of node 10 (cold, 290 its), the
# polish variants only (mode 34: HG and HGR, lean), with 32 exact passes and
# 48 reduce passes + 1 exact, evaluated at every check within 20 of the
# cutoff (no other checks). Compare with E1 (8 exact, 12 + 1 reduce).
# Usage: run_e5.sh <core>
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
f=$D/logs/e5_c37_passes.log
grep -q '^RESULT' "$f" 2>/dev/null && exit 0
DNN_CERTX=34 DNN_CERTX_PASSES=32 DNN_CERTX_RPASSES=48 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=20 \
  taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 \
  -x "39:27,30:48,25:49,37:45" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
