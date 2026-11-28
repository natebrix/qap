#!/bin/bash
# run.sh PROBLEM UB ARM: PROBLEM to the end under nug30's recipe with strong branching at depth 1 (sb.early = 2) and
# sb.reuse = ARM (0: off; 1: fathoms and bounds; 2: and the children's states), in one process.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); P=$1; UB=$2; R=$3; T=$P.r$R
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$D/../../../latest"
rm -f "$D/runs/$T.dnn.txt" "$D/runs/$T.sb.txt"
QAP_DNN_TRACE="$D/runs/$T.dnn.txt" QAP_SB_TRACE="$D/runs/$T.sb.txt" "${BIN:-./bnbdrv}" -p "$D/par/r$R.par" -u "$UB" -o "$D/runs/$T.state" ../qaplib/$P.dat > "$D/runs/$T.out.txt" 2>&1
echo "$T exit $? $(date -u +%T)" >> "$D/runs/done.txt"
