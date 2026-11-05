#!/bin/bash
# The node-10 subtree through the coordinator, two arms at the same time on
# one machine (README.md). Usage, from anywhere in the repository:
#
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh base zone
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh warmy certx
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh zone2 base2   (replicates)
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh full ship2    (all levers; ship replicate)
#
# Arms (each on top of base unless stated):
#   base   the production configuration of coord_search/run_coord_search.sh:
#          rule 9, bundle D, parent-Y warm start, pool inheritance
#   zone   base + QAP_DNN_CUT_ZONE=0.005: misses up to 244 short (node 10
#          among them) run a cut phase, so their children inherit a pool
#   warmy  zone + Y kept between cut rounds (QAP_DNN_CUT_WARMY=1), 12 rounds
#          of 200 iterations, quiet round at most 1,000, first Polyak step
#          damped for phases that start from an inherited pool
#   certx  zone + the certificate polish within 10 units of the cutoff
#   ship   zone + the parent's DNN state carried across units in the state
#          files (QAP_DNN_SHIP=$SHIP_BITS, default 32), so no node starts cold
#          because its parent was bounded in another unit
#   full   zone + warmy's phase settings + ship: every lever that passed
#   zone2, base2, ship2  replicates of zone, base and ship: the run-to-run noise of one
#          configuration on one machine (unit splits, steals, marginal nodes)
#
# Each arm gets SLOTS solver processes (default: half the performance cores,
# at least 2). An arm that has finished is skipped, and an interrupted arm
# resumes where it stopped: run the same command again. On macOS the machine
# is kept awake while the script runs (caffeinate; keep the lid open or the
# power and an external display connected).
#
# Results: runs/ARM/ (trace.txt, results.txt, coordinator.log, worker.log,
# hosts.txt, env.txt), and summary.txt comparing every finished arm with base.
# The stores and work directories (runs/ARM/backend, work, bin) are ignored by
# git.
set -u
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../.." && pwd)
RUNS=${PAB_RUNS:-$D/runs}

[ $# -ge 1 ] || { sed -n 2,30p "$0"; exit 2; }
for arm in "$@"; do
  case $arm in base|zone|warmy|certx|ship|full|zone2|base2|ship2) ;; *) echo "unknown arm '$arm' (base, zone, warmy, certx, ship, full, zone2, base2, ship2)"; exit 2;; esac
done

if [ "$(uname -s)" = Darwin ]; then
  P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
  CPU="$(sysctl -n machdep.cpu.brand_string) P-cores $P"
  command -v caffeinate > /dev/null && { caffeinate -ims -w $$ & }
else
  P=$(nproc)
  CPU="$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') amx_cpus $(grep -c amx_bf16 /proc/cpuinfo)"
fi
SLOTS=${SLOTS:-$(( P / 2 > 2 ? P / 2 : 2 ))}

cd "$R" || exit 1
[ -x latest/bnbdrv ] || { echo "latest/bnbdrv is missing: run 'make solver' first"; exit 1; }
PY=${PYTHON:-python3}

# the production settings (coord_search/run_coord_search.sh)
INSTANCE=${PAB_INSTANCE:-qaplib/wil50.dat}
UB=${PAB_UB:-48816}
EXTRA=${PAB_EXTRA:--F 39:27,30:48,25:49 -x 0:4:1:300:1}
GRID=${PAB_GRID:-5:10}
BUDGET=${PAB_BUDGET:-3600}
UNIT=${PAB_UNIT:-1800}

arm_env() {
  export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
  export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=$GRID QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
  export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
  export QAP_DNN_WARM_Y=2 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
  case $1 in
    zone|zone2) export QAP_DNN_CUT_ZONE=0.005 ;;
    warmy) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 ;;
    full)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} ;;
    ship|ship2) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_SHIP=${SHIP_BITS:-32} ;;
    certx) export QAP_DNN_CUT_ZONE=0.005 DNN_CERTX=32 DNN_CERTX_USE=1 DNN_CERTX_QUIET=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=10 ;;
  esac
}

run_arm() {
  local arm=$1 A=$RUNS/$1
  local B=local:$A/backend
  mkdir -p "$A"
  if [ -f "$A/backend/store/runs/$arm/results.txt" ]; then
    echo "$arm: finished earlier"
    return 0
  fi
  arm_env "$arm"
  if [ ! -x "$A/bin/bnbdrv" ]; then
    # a private copy, so that a rebuild during the run does not change it
    if { [ "$arm" = warmy ] || [ "$arm" = full ]; } && ! grep -aq QAP_DNN_CUT_DAMP1 latest/bnbdrv; then
      echo "warmy: this bnbdrv has no first-step damping yet: git pull, make solver, run again"; return 1
    fi
    if { [ "$arm" = ship ] || [ "$arm" = ship2 ] || [ "$arm" = full ]; } && ! grep -aq QAP_DNN_SHIP latest/bnbdrv; then
      echo "ship: this bnbdrv cannot carry the DNN state yet: git pull, make solver, run again"; return 1
    fi
    mkdir -p "$A/bin" && cp latest/bnbdrv "$A/bin/bnbdrv"
    echo "$(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)')" > "$A/bin/commit.txt"
  fi
  env | grep -E '^(QAP_|DNN_)' | sort > "$A/env.txt"
  echo "$(date -u +%FT%TZ) $CPU slots $SLOTS $(cat "$A/bin/commit.txt")" >> "$A/hosts.txt"
  if [ ! -f "$A/backend/ledger.db" ]; then
    $PY -m coord.qapctl init "$arm" --backend "$B" --instance "$INSTANCE" --params "$R/latest/PARAM/nug_deep_xbar9_all.par" \
      --ub "$UB" --slots "$SLOTS" --extra "$EXTRA" --trace --solver "$A/bin/bnbdrv" \
      --budget "$BUDGET" --target-unit-cpu "$UNIT" --steal-after 10 --lease 600 --heartbeat 60 --ckpt-upload 300 \
      --checkpoint-every 60 --status-every 60 >> "$A/coordinator.log" 2>&1 || { echo "$arm: init failed (see $A/coordinator.log)"; return 1; }
  else
    $PY -m coord.qapctl requeue-leased "$arm" --backend "$B" >> "$A/coordinator.log" 2>&1
  fi
  $PY -m coord.qapctl worker "$arm" --backend "$B" --slots "$SLOTS" --solver "$A/bin/bnbdrv" \
    --workdir "$A/work" >> "$A/worker.log" 2>&1 &
  echo $! > "$A/worker.pid"
  $PY -m coord.qapctl coordinator "$arm" --backend "$B" >> "$A/coordinator.log" 2>&1 &
  echo $! > "$A/coordinator.pid"
  wait "$(cat "$A/coordinator.pid")"
  local rc=$?
  wait "$(cat "$A/worker.pid")"
  rm -f "$A/worker.pid" "$A/coordinator.pid"
  if [ -f "$A/backend/store/runs/$arm/results.txt" ]; then
    find "$A/backend/store/runs/$arm/units" -name trace.txt -exec cat {} + > "$A/trace.txt"
    cp "$A/backend/store/runs/$arm/results.txt" "$A/results.txt"
    echo "$arm: done, $(tail -1 "$A/results.txt")"
  else
    echo "$arm: stopped before the end (coordinator exit $rc); run the same command again to resume"
  fi
}

echo "$(date '+%F %T') arms: $* | $CPU | $SLOTS slots per arm | logs in $RUNS/ARM/ | stop with: kill $$"
mkdir -p "$RUNS" && echo $$ > "$RUNS/pair.pid"
# kill PID (TERM) or ^C: the workers drain (open nodes go back to their queues
# with their checkpoints) and the coordinators stop; the same command resumes
stop() {
  echo "$(date '+%F %T') stopping: the arms drain and can be resumed with the same command"
  for f in "$RUNS"/*/worker.pid; do [ -f "$f" ] && kill -TERM "$(cat "$f")" 2>/dev/null; done
  sleep 5
  for f in "$RUNS"/*/coordinator.pid; do [ -f "$f" ] && kill -TERM "$(cat "$f")" 2>/dev/null; done
}
trap stop INT TERM
pids=()
for arm in "$@"; do
  ( trap - INT TERM; run_arm "$arm" ) &
  pids+=($!)
done
for p in "${pids[@]}"; do
  while kill -0 "$p" 2>/dev/null; do wait "$p"; done   # wait again after a trapped signal
done
rm -f "$RUNS/pair.pid"
# the summary goes next to the runs it describes (summary.txt for runs/)
SUMMARY=$D/summary.txt; [ "$RUNS" != "$D/runs" ] && SUMMARY=$RUNS/summary.txt
$PY "$D/summarize.py" "$RUNS" | tee "$SUMMARY"
echo
echo "To share the results: git add experiments/wil50_feasibility/pipeline_ab && git commit -m 'pipeline A/B: $*' && git push"
