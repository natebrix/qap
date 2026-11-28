#!/bin/bash
# s42/run.sh NODE [LINES]: one DNN strong-branching firing at a depth-2 node of sko42 (NODE = its two fixes, e.g.
# 22:0,24:41) under the run of record's recipe (sko42-261124, infra/phase1/recipes/sko42_d2r10_sb.par) with
# sb.depths = 2:2, every child of every candidate solved, and LINES (default c35,c6, the lines rule 10 took at
# that run's heavy depth-2 nodes) added to the candidates. Cold start (-F), stops at depth 2 (-d 2). Resumable.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); node=$1; lines=${2:-c35,c6}; name=${node//:/-}; name=${name//,/_}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$D/../../../latest" && QAP_SB_LINES="$lines" QAP_SB_TRACE="$D/runs/$name.sb.txt" QAP_SB_PROGRESS="$D/runs/$name.progress.txt" \
  QAP_DNN_TRACE="$D/runs/$name.dnn.txt" "${BIN:-./bnbdrv}" -p "$D/par/d2.par" -s SYM/sko42t.sym -u 15814 -F "$node" -d 2 \
  -o "$D/runs/$name.state" ../qaplib/sko42t.dat >> "$D/runs/$name.out.txt" 2>&1
echo "$name exit $? $(date -u +%FT%TZ)" >> "$D/runs/done.txt"
