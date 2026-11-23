#!/bin/bash
# ucheap/sko42root: the root's DNN strong-branching firing on sko42 under ../suggest/sko42_d2r10_sb.par
# (every row, the 8th-power score), then stop (-d 0): which row does it take, and what does it cost?
# sko42-261120b branched the root on r36 (rule 9 without marginals: rule 4's lookahead).
# Resumable: QAP_SB_PROGRESS keeps every child solve.   run.sh   (results in runs/)
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$D/runs"
cd "$BUILD/latest" && QAP_SB_TRACE="$D/runs/sb.txt" QAP_SB_PROGRESS="$D/runs/progress.txt" QAP_DNN_TRACE="$D/runs/dnn.txt" \
  ./bnbdrv -p "$D/../suggest/sko42_d2r10_sb.par" -s SYM/sko42t.sym -u 15814 -d 0 -o "$D/runs/out.state" \
  ../qaplib/sko42t.dat >> "$D/runs/out.txt" 2>&1
echo "exit $? $(date -u +%T)" >> "$D/runs/done.txt"
