#!/bin/bash
# d2/run.sh NODE [LINES]: one DNN strong-branching firing at a depth-2 node of wil50 (NODE = its two fixes, e.g.
# 21:20,45:0) under the live run's recipe with sb.depths = 2:2, every child of every candidate solved (sb.early
# off), and LINES (default r31, the line wil50-261126 took at all 162 of its depth-2 nodes) added to the
# candidates. The node starts cold here (-F), where the run's started from its parent's state and pool.
# Resumable (QAP_SB_PROGRESS). Output in runs/NAME.*
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); node=$1; lines=${2:-r31}; name=${node//:/-}; name=${name//,/_}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$D/../../../latest" && QAP_SB_LINES="$lines" QAP_SB_TRACE="$D/runs/$name.sb.txt" QAP_SB_PROGRESS="$D/runs/$name.progress.txt" \
  QAP_DNN_TRACE="$D/runs/$name.dnn.txt" "${BIN:-./bnbdrv}" -p "$D/par/d2.par" -s SYM/wil50t.sym -u 48818 -F "$node" -d 2 \
  -o "$D/runs/$name.state" ../qaplib/wil50t.dat >> "$D/runs/$name.out.txt" 2>&1
echo "$name exit $? $(date -u +%FT%TZ)" >> "$D/runs/done.txt"
