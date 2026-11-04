#!/bin/bash
# E2: the certificate lag in the cut phase. One node through the production
# phase protocol (an uncut round 0 of <first> iterations, the search's own
# give-up point for that node; rounds of 300 iterations with 1,000 new
# regional triangles each, Polyak target best + 50 uncapped, stall 3, then a
# 1,500-iteration quiet round from the best round's multipliers), with the
# certificate variants evaluated every <every> checks of every solve (no
# window). The questions: how much of the end-of-round lag (certified value
# against the quiet round's value) do the variants recover at the same
# iterate, and would they fathom the node earlier (a round end instead of the
# quiet round)? With bit 16 in the mode, M of the last evaluated check is
# dumped to <dumpdir>/<name> for scripts/child_gl.py (question 3).
# dnncuts has no primal gate on the quiet round (the search runs it only when
# the best round's primal is above the cutoff).
# Usage: run_e2.sh <core> <mode> <name> <fixes> <first> <every> [rounds] [dumpdir]
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; mode=$2; name=$3; fixes=$4; first=$5; every=$6; rounds=${7:-12}; dd=${8:-/tmp}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
f=$D/logs/e2_${name}_m${mode}.log
grep -q '^FINAL' "$f" 2>/dev/null && exit 0
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core" > "$f.tmp"
DNN_CERTX=$mode DNN_CERTX_EVERY=$every DNN_CERTX_WINDOW=-1e9 DNN_CERTX_DUMP=$dd/$name \
  taskset -c "$core" "$R/experiments/dnn_cuts/dnncuts" -x "$fixes" -grid 5:10 -u 48816 -T 1 -d 3 \
  -I "$first" -i 300 -r "$rounds" -k 1000 -delta 50 -nocap -stall 3 -freeze 0:1500 -v \
  "$R/qaplib/wil50.dat" >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
