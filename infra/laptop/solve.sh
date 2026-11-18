#!/bin/bash
# A whole solve on one machine (the laptop) through the coordinator's local mode, with the
# settings start.sh would launch on AWS (infra/phase1/README.md, Starting a run), and resumable:
#
#   caffeinate -i infra/laptop/solve.sh PROBLEM
#
# The defaults, as start.sh's: the params file recipes/defaults.txt names for PROBLEM (else
# recipes/PROBLEM.par, else experiments/estimates/PROBLEM/PROBLEM.par); the data file
# qaplib/PROBLEMt.dat with latest/SYM/PROBLEMt.sym where the grid's symmetry file exists, else
# qaplib/PROBLEM.dat; the incumbent the best known (QAPSOLN/) + the step of the objective's grid
# (objstep.py); the unit sizing of recipes/PROBLEM.sh if there is one. The run is named PROBLEM-laptop
# and kept in OUT (default infra/laptop/runs/PROBLEM/, git-ignored but for its results).
#
# Interrupted (the lid closed, Ctrl-C, a reboot): run the same command again. The run goes on from
# its store: the units a dead worker held go back to the queue (qapctl requeue-leased), and the
# coordinator replays what was merged. Only a finished run is not repeated.
#
# Environment: SLOTS (default: the cores less 2), OUT, MIN_FREE_GB (10: the run is stopped below it),
# and for tests PARAMS (a params file relative to the repository), UB, COORD (unit sizing).
#
# ESTIMATE=SOLVE_S:DIVE_S[:B] makes it an estimate run (docs/estimates_plan.md; coord/dives.py): the
# solve for SOLVE_S seconds of its wall time, then paused, and DIVE_S seconds of survivor dives per
# open unit (weights exp(B x), default 0.4) as run PROBLEM-est-dives; OUT/estimate.txt at the end.
# The run is PROBLEM-est (OUT default infra/laptop/runs/PROBLEM-est/), and the same command resumes
# it at whatever step it was. The solve stays paused in the store: `qapctl resume` goes on with it.
set -uo pipefail
cd "$(dirname "$0")/../.."
R=$(pwd)
P=${1:-}
[ -n "$P" ] || { awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 2; }
[[ $P =~ ^[a-z]+[0-9]+[a-z]?$ ]] && [ -f "qaplib/$P.dat" ] || { echo "no QAPLIB problem $P (qaplib/$P.dat)"; exit 2; }
PY=${PYTHON:-python3}
SLOTS=${SLOTS:-$(( $( (sysctl -n hw.ncpu 2> /dev/null || nproc) ) - 2 ))}
[ "$SLOTS" -ge 1 ] || SLOTS=1
EST=${ESTIMATE:-}
EST_SOLVE= EST_DIVE= EST_B=
if [ -n "$EST" ]; then
  IFS=: read -r EST_SOLVE EST_DIVE EST_B <<< "$EST"
  EST_B=${EST_B:-0.4}
  [ -n "$EST_SOLVE" ] && [ -n "$EST_DIVE" ] || { echo "ESTIMATE must be SOLVE_S:DIVE_S[:B]"; exit 2; }
  OUT=${OUT:-$R/infra/laptop/runs/$P-est}
  RUN=$P-est
else
  OUT=${OUT:-$R/infra/laptop/runs/$P}
  RUN=$P-laptop
fi
MIN_FREE_GB=${MIN_FREE_GB:-10}
ST=$OUT/store
BK="local:$ST"
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$OUT"
say() { echo "$(date '+%m-%d %H:%M:%S') $*" | tee -a "$OUT/run.log"; }
free_gb() { df -Pk "$OUT" | awk 'NR == 2 {printf "%d", $4 / 1048576}'; }

# --- the problem's defaults (as start.sh)
BEST=""; [ -f "QAPSOLN/$P.sln" ] && BEST=$(awk 'NR == 1 {print $2; exit}' "QAPSOLN/$P.sln" | tr -d '\r[:space:]')
[[ $BEST =~ ^[0-9]+$ ]] || BEST=""
DAT=qaplib/$P.dat; SYM=""
if [ -f "qaplib/${P}t.dat" ] && [ -f "latest/SYM/${P}t.sym" ]; then DAT=qaplib/${P}t.dat; SYM=SYM/${P}t.sym
elif [ -f "latest/SYM/$P.sym" ]; then SYM=SYM/$P.sym; fi
PAR=${PARAMS:-$(awk -v p="$P" '!/^#/ && $1 == p {print $2; exit}' infra/phase1/recipes/defaults.txt 2> /dev/null)}
if [ -z "$PAR" ] || [ ! -f "$PAR" ]; then
  PAR=""
  for c in "infra/phase1/recipes/$P.par" "experiments/estimates/$P/$P.par" "experiments/${P}_feasibility/fullrun/fullrun.par"; do
    [ -f "$c" ] && { PAR=$c; break; }
  done
fi
[ -n "$PAR" ] || { echo "no params file for $P: give one with PARAMS=path (relative to the repository)"; exit 2; }
STEP=$($PY infra/phase1/objstep.py "$DAT" | awk '{print $1}') || exit 1
UBV=${UB:-}; [ -z "$UBV" ] && [ -n "$BEST" ] && UBV=$((BEST + STEP))
RCOORD=""
[ -f "infra/phase1/recipes/$P.sh" ] && RCOORD=$(COORD=; source "infra/phase1/recipes/$P.sh" > /dev/null; echo "$COORD")
COORDV=${COORD:-${RCOORD:-"--budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300 --trace"}}

res=$OUT/results.txt
if [ -s "$res" ]; then say "$P: finished before ($res)"; exit 0; fi
if [ -n "$EST" ] && [ -s "$OUT/estimate.txt" ]; then say "$P: estimated before ($OUT/estimate.txt)"; exit 0; fi
say "$P: $DAT${SYM:+, -s $SYM}, params $PAR, incumbent ${UBV:-none} (best known ${BEST:-?} + $STEP), $SLOTS slots, $(free_gb) GB free, commit $(git rev-parse --short HEAD)"
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { say "build failed: $OUT/build.log"; exit 1; }

# --- stop the run when the disk runs low (parents' states and unit attempts live in OUT)
( while sleep 120; do
    if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
      say "disk below $MIN_FREE_GB GB: stopping the run (run this again once there is room)"
      pkill -TERM -f "coord.qapctl (local|coordinator|worker|estimate-flow) $RUN( |$)" 2> /dev/null
      break
    fi
  done ) &
watch=$!
trap 'kill $watch 2> /dev/null' EXIT

status=$($PY infra/phase1/bk.py "$BK" status "$RUN" 2> /dev/null || echo none)
if [ -n "$EST" ]; then
  # an estimate run: init (a new run), then a worker that follows the run into its dives run and
  # qapctl estimate-flow in the coordinator's place; both pick up where a stopped attempt left off
  if [ "$status" = none ]; then
    say "$P: a new estimate run (the root first, in qapctl init; status: $OUT/coord.log)"
    $PY -m coord.qapctl init "$RUN" --backend "$BK" --instance "$R/$DAT" --params "$R/$PAR" \
      ${UBV:+--ub "$UBV"} --slots "$SLOTS" ${SYM:+--extra "-s $R/latest/$SYM"} $COORDV \
      --visibility private >> "$OUT/coord.log" 2>&1 || { say "init failed: $OUT/coord.log"; exit 1; }
  else
    say "$P: resuming the estimate run ($status, dives run $($PY infra/phase1/bk.py "$BK" status "$RUN-dives"))"
    $PY -m coord.qapctl requeue-leased "$RUN" --backend "$BK" >> "$OUT/coord.log" 2>&1
    [ "$($PY infra/phase1/bk.py "$BK" status "$RUN-dives")" = none ] \
      || $PY -m coord.qapctl requeue-leased "$RUN-dives" --backend "$BK" >> "$OUT/coord.log" 2>&1
  fi
  $PY -m coord.qapctl worker "$RUN" --backend "$BK" --slots "$SLOTS" --workdir "$OUT/work" --follow \
    >> "$OUT/worker.log" 2>&1 &
  wpid=$!
  $PY -m coord.qapctl estimate-flow "$RUN" --backend "$BK" --solve-wall "$EST_SOLVE" --dive-budget "$EST_DIVE" \
    --dive-b "$EST_B" >> "$OUT/coord.log" 2>&1
  rc=$?
  wait $wpid 2> /dev/null
  if [ $rc = 0 ] && [ -s "$ST/store/runs/$RUN-dives/results.txt" ]; then
    cp "$ST/store/runs/$RUN-dives/results.txt" "$OUT/estimate.txt"
    cp "$ST/store/runs/$RUN-dives/estimate.json" "$OUT/estimate.json" 2> /dev/null
    say "$P: $(grep -E '^estimated tree' "$OUT/estimate.txt")"
    say "$P: the estimate is in $OUT/estimate.txt; the solve $RUN is paused in the store"
    exit 0
  elif [ $rc = 0 ] && [ -s "$ST/store/runs/$RUN/results.txt" ]; then
    cp "$ST/store/runs/$RUN/results.txt" "$res"
    say "$P: the tree closed inside the solve's budget: $res"
    exit 0
  fi
  say "$P: the estimate did not finish (exit $rc): run the same command again to go on; logs in $OUT"
  exit 1
fi
case $status in
  none)
    say "$P: a new run (the root first, in qapctl init; status: $OUT/coord.log)"
    $PY -m coord.qapctl local "$RUN" --backend "$BK" --instance "$R/$DAT" --params "$R/$PAR" \
      ${UBV:+--ub "$UBV"} --slots "$SLOTS" ${SYM:+--extra "-s $R/latest/$SYM"} $COORDV \
      --visibility private --workdir "$OUT/work" >> "$OUT/coord.log" 2>&1
    rc=$? ;;
  finished) rc=0 ;;
  stopped) say "$P: the run was stopped (qapctl stop): not resumed"; exit 1 ;;
  *)
    say "$P: resuming ($status): the units a dead worker held go back to the queue"
    $PY -m coord.qapctl requeue-leased "$RUN" --backend "$BK" >> "$OUT/coord.log" 2>&1
    $PY -m coord.qapctl worker "$RUN" --backend "$BK" --slots "$SLOTS" --workdir "$OUT/work" >> "$OUT/worker.log" 2>&1 &
    wpid=$!
    $PY -m coord.qapctl coordinator "$RUN" --backend "$BK" >> "$OUT/coord.log" 2>&1
    rc=$?
    kill -TERM $wpid 2> /dev/null; wait $wpid 2> /dev/null ;;
esac

if [ $rc = 0 ] && [ -s "$ST/store/runs/$RUN/results.txt" ]; then
  cp "$ST/store/runs/$RUN/results.txt" "$res"
  case " $COORDV " in *" --trace "*)   # the node and DNN traces of a --trace run
    $PY -m coord.qapctl export "$RUN" --backend "$BK" --out "$OUT/export" >> "$OUT/coord.log" 2>&1 \
      || say "export failed: the store is kept in $ST" ;;
  esac
  say "$P: $(grep -E 'Nodes:|coordinator:' "$res" | tr -s ' ' | tr '\n' ' ' | cut -c1-240)"
  got=$(grep -oE 'best objective [0-9]+' "$res" | tail -1 | awk '{print $3}')
  found=$(grep -c "New objective value: ${BEST:-x} " "$OUT/coord.log" 2> /dev/null || true)
  if [ -n "$BEST" ] && [ "$got" = "$BEST" ]; then
    say "$P: the tree closed at the best known value $BEST$([ "${found:-0}" -gt 0 ] && echo ", found by the search and verified by the coordinator")"
  else
    say "$P: ended at ${got:-?}, NOT the best known ${BEST:-?}: read $res before anything else"
  fi
else
  say "$P: not finished (exit $rc): run the same command again to go on; logs in $OUT"
  exit 1
fi
