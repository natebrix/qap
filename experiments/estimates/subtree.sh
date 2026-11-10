#!/bin/bash
# The complete subtree of one node, through the coordinator, with an instance's
# current settings (README.md here):
#
#   experiments/estimates/subtree.sh SETTINGS NAME SPEC [SLOTS]
#   experiments/estimates/subtree.sh experiments/wil50_feasibility/fullrun/settings.sh d3_220 21:22,31:23,40:15
#
# SETTINGS is a settings file (the full wil50 run's, or INSTANCE/settings.sh here);
# its params file, instance, symmetry file and incumbent are used as they are, after
# its fullrun_check. SPEC is the node's assignment in that instance's indices (-F).
# SLOTS defaults to the performance cores (all of them: one subtree at a time).
# Output in RUNDIR/NAME (RUNDIR: the settings file's directory's subtrees/, or
# $EST_RUNS): trace.txt (one line per DNN-bounded node), results.txt, summary.txt.
# ^C or kill drains (open nodes go back to the queue with their checkpoints); the
# same command resumes.
set -u
[ $# -ge 3 ] || { sed -n '2,/^set -u$/p' "$0" | sed '$d'; exit 2; }
SETF=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); NAME=$2; SPEC=$3
R=$(cd "$(dirname "$0")/../.." && pwd)
cd "$R" || exit 1
. "$SETF" || exit 1
fullrun_check || exit 1
if [ "$(uname -s)" = Darwin ]; then
  P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
  CPU="$(sysctl -n machdep.cpu.brand_string) P-cores $P"
  command -v caffeinate > /dev/null && { caffeinate -ims -w $$ & }
else
  P=$(nproc); CPU="$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') amx_cpus $(grep -c amx_bf16 /proc/cpuinfo)"
fi
SLOTS=${4:-$P}
RUNS=${EST_RUNS:-$(dirname "$SETF")/subtrees}
A=$RUNS/$NAME; B=local:$A/backend; RUN=$NAME
mkdir -p "$A/bin"
[ -x "$A/bin/bnbdrv" ] || { cp latest/bnbdrv "$A/bin/bnbdrv"; echo "$(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)')" > "$A/bin/commit.txt"; }
EXTRA="-F $SPEC${FULLRUN_SYM:+ -s $FULLRUN_SYM}"
{ echo "settings $SETF"; echo "spec $SPEC"; echo "extra $EXTRA"; env | grep -E '^(QAP_|DNN_)' | sort; } > "$A/env.txt"
echo "$(date -u +%FT%TZ) $CPU slots $SLOTS $(cat "$A/bin/commit.txt")" >> "$A/hosts.txt"
PY=${PYTHON:-python3}
if [ ! -f "$A/backend/ledger.db" ]; then
  $PY -m coord.qapctl init "$RUN" --backend "$B" --instance "$FULLRUN_INSTANCE" --params "$FULLRUN_PAR" \
    --ub "$FULLRUN_UB" --slots "$SLOTS" --extra "$EXTRA" --trace --solver "$A/bin/bnbdrv" \
    --budget 3600 --target-unit-cpu 1800 --steal-after 10 --lease 600 --heartbeat 60 --ckpt-upload 300 \
    --checkpoint-every 60 --status-every 60 >> "$A/coordinator.log" 2>&1 || { echo "init failed (see $A/coordinator.log)"; exit 1; }
else
  $PY -m coord.qapctl requeue-leased "$RUN" --backend "$B" >> "$A/coordinator.log" 2>&1
fi
echo "$(date '+%F %T') $NAME: $SPEC with $(basename "$(dirname "$SETF")")'s settings, $SLOTS slots, logs in $A | stop with: kill $$"
$PY -m coord.qapctl worker "$RUN" --backend "$B" --slots "$SLOTS" --solver "$A/bin/bnbdrv" --workdir "$A/work" >> "$A/worker.log" 2>&1 &
WPID=$!
$PY -m coord.qapctl coordinator "$RUN" --backend "$B" >> "$A/coordinator.log" 2>&1 &
CPID=$!
stop() { echo "$(date '+%F %T') stopping: draining; the same command resumes"; kill -TERM $WPID 2>/dev/null; sleep 5; kill -TERM $CPID 2>/dev/null; }
trap stop INT TERM
while kill -0 $CPID 2>/dev/null; do wait $CPID; done
while kill -0 $WPID 2>/dev/null; do wait $WPID; done
if [ -f "$A/backend/store/runs/$RUN/results.txt" ]; then
  find "$A/backend/store/runs/$RUN/units" -name trace.txt -exec cat {} + > "$A/trace.txt"
  cp "$A/backend/store/runs/$RUN/results.txt" "$A/results.txt"
  $PY "$R/experiments/estimates/subtree_summary.py" "$A" | tee "$A/summary.txt"
else
  echo "$NAME: stopped before the end; run the same command again to resume"
fi
