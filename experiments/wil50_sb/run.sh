#!/bin/bash
# wil50_sb/run.sh [NODE]: one DNN strong-branching firing on wil50 under par/sb.par (the estimate's recipe,
# wil50_d2r10_g150xl2.par, with strong branching at the root and depth 1), then stop: at the root (no NODE,
# -d 0) or at a depth-1 node (-F NODE -d 1). Which line does it take, by what margin, at what cost?
# Resumable: QAP_SB_PROGRESS keeps every child solve. Output in runs/NAME.*
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); BUILD=${BUILD:-/home/user/wt-rule10}
node=${1:-}; name=${node:-root}; name=${name//:/-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$BUILD/latest" && QAP_SB_TRACE="$D/runs/$name.sb.txt" QAP_SB_PROGRESS="$D/runs/$name.progress.txt" \
  QAP_DNN_TRACE="$D/runs/$name.dnn.txt" QAP_NODE_TRACE="$D/runs/$name.nodes.txt" \
  ./bnbdrv -p "$D/par/sb.par" -s SYM/wil50t.sym -u 48818 ${node:+-F "$node"} -d $([ -n "$node" ] && echo 1 || echo 0) \
  -o "$D/runs/$name.state" ../qaplib/wil50t.dat >> "$D/runs/$name.out.txt" 2>&1
echo "$name exit $? $(date -u +%FT%TZ)" >> "$D/runs/done.txt"
