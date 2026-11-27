#!/bin/bash
# nug30 record attempt on the laptop (Nov 19): g5e (the Nov 16 record, bound_lab/nug30_phase) against
# one change at a time. Each arm is the whole search from incumbent 6124 with traces on, as the record
# was run (bound_lab/laptop/nug30.sh); it must end at 6124 with no open node.
#   experiments/nug30_record/run.sh ARM[@JOBS] ...    e.g.  g5e g5e@9 g5e@10 b3r40 two gu83 wcg10
# ARM is par/ARM.par; JOBS defaults to 8. A repeat of a label gets a suffix (g5e, g5e-2, ...).
# Each run's directory has the log, the traces, the commit and the load before it started.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../.."
R=$(pwd)
PY=${PYTHON:-python3}
OUT=${OUT:-$H/runs}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "$OUT"
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { echo "build failed: $OUT/build.log"; exit 1; }
for spec in "$@"; do
  arm=${spec%@*}; jobs=8; [ "$spec" != "$arm" ] && jobs=${spec#*@}
  label=$arm; [ "$jobs" != 8 ] && label=$arm-j$jobs
  D=$OUT/$label; n=2; while [ -e "$D" ]; do D=$OUT/$label-$n; n=$((n + 1)); done
  mkdir -p "$D"
  W=$(mktemp -d "${TMPDIR:-/tmp}/n30-$label-XXXXXX")
  { git rev-parse --short HEAD; git status --porcelain latest/ | wc -l | sed 's/^ */latest\/ uncommitted files: /'; uptime; } > "$D/env.txt"
  echo "$(date +%H:%M:%S) $(basename "$D"): $arm.par, $jobs jobs"
  QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt \
    "$PY" bench/parallel.py qaplib/nug30.dat --jobs "$jobs" --ub 6124 \
      --params "$H/par/$arm.par" --workdir "$W" \
      --checkpoint "$D/final.state" --checkpoint-every 100000 > "$D/log.txt" 2> "$D/err.txt" \
    || echo "$(basename "$D"): exit $? (see $D/err.txt)"
  rm -rf "$W" "$D/final.state"
  grep 'parallel: .* workers' "$D/log.txt" | grep -q 'best objective 6124' || echo "$(basename "$D"): did not end at 6124"
  grep -h 'parallel: .* workers' "$D/log.txt" | tail -1
done
