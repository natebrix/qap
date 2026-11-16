#!/bin/bash
# nug30 on the laptop: the record recipe against today's best arm (bound_lab/nug30_phase).
#   base   the record recipe (rule 9, DNN at depths 1-7): 3.0 min on the M2 Pro on Nov 8
#   g5e    long phases at depth 1, short below, the phase's own penalty, early entry:
#          1.82x less CPU and wall than base on the cloud box (needs this branch's solver: cut.giveup_depth)
# From the repository's root:   experiments/bound_lab/laptop/nug30.sh
# Options (environment): ARMS ("base g5e"), JOBS (8), OUT (runs_nug30/ here). Each arm is the whole
# search from incumbent 6124; it must end with no open node.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../.."
R=$(pwd)
ARMS=${ARMS:-"base g5e"}
JOBS=${JOBS:-8}
OUT=${OUT:-$H/runs_nug30}
PY=${PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "$OUT"
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { echo "build failed: $OUT/build.log"; exit 1; }
for arm in $ARMS; do
  D=$OUT/$arm; rm -rf "$D"; mkdir -p "$D"
  W=$(mktemp -d "${TMPDIR:-/tmp}/n30-$arm-XXXXXX")
  echo "$(date +%H:%M:%S) $arm: running, $JOBS jobs"
  QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt \
    "$PY" bench/parallel.py qaplib/nug30.dat --jobs "$JOBS" --ub 6124 \
      --params "$R/experiments/bound_lab/nug30_phase/par/$arm.par" --workdir "$W" \
      --checkpoint "$D/final.state" --checkpoint-every 100000 > "$D/log.txt" 2> "$D/err.txt" \
    || echo "$arm: exit $? (see $D/err.txt)"
  rm -rf "$W" "$D/final.state"
done
"$PY" "$H/nug30_summary.py" "$OUT" $ARMS | tee "$OUT/summary.txt"
