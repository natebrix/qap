#!/bin/bash
# wil50_night: one complete wil50 subtree under one arm, cold from -F, one process (README.md).
#   run.sh NAME SPEC ARM [ENV=VAL ...]     e.g. run.sh a220 21:22,31:23,40:15 base
# par/ARM.par; results in runs/NAME_ARM/ (log.txt, dnn.txt, nodes.txt). Incumbent 48818, the next run's.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
name=$1; spec=$2; arm=$3; shift 3
O=$D/runs/${name}_$arm; mkdir -p "$O"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for kv in "$@"; do export "$kv"; done
echo "$(date -u +%FT%TZ) $spec $arm $* build $(git -C "$BUILD" rev-parse --short HEAD)" > "$O/cmd.txt"
cd "$BUILD/latest" && QAP_DNN_TRACE="$O/dnn.txt" QAP_DNN_TRACE_HIST=1 QAP_NODE_TRACE="$O/nodes.txt" \
  ./bnbdrv -p "$D/par/$arm.par" -s SYM/wil50t.sym -u 48818 -F "$spec" -o "$O/out.state" ../qaplib/wil50t.dat > "$O/log.txt" 2> "$O/err.txt"
echo "exit $? $(date -u +%T)" >> "$O/log.txt"
echo "$(date -u +%T) ${name}_$arm: $(grep -E 'Nodes:|Total Time' "$O/log.txt" | tr -s ' ' | tr '\n' ' ')" >> "$D/runs.log"
