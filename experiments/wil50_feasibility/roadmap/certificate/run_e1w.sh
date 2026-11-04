#!/bin/bash
# E1w: as E1 but the children start from node 10's final state (Z and the
# conditioned Y, QAP_DNN_WARM_Y=2 = dnntest -Q b), as the coordinator's
# in-memory warm starts do. Node 10 is solved once (300 iterations, give-up),
# not instrumented (DNN_CERTX_CALLS skips call 0).
# Usage: run_e1w.sh <core> <mode> <child> [<child> ...]   (child = X for X:45)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; mode=$2; shift 2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
kids=""; for x in "$@"; do kids="$kids${kids:+;}$x:45"; done
f=$D/logs/e1w_m${mode}_$(echo "$@" | tr ' ' '_').log
grep -q '^RESULT node=1' "$f" 2>/dev/null && exit 0
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core" > "$f.tmp"
DNN_CERTX=$mode DNN_CERTX_CALLS=$(( (1 << 40) - 2 )) DNN_CERTX_EVERY=${EVERY:-5} DNN_CERTX_WINDOW=${WIN:-50} \
  taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 \
  -x "39:27,30:48,25:49" -K "$kids" -Q b "$R/qaplib/wil50.dat" >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
