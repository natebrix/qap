#!/bin/bash
# The node-10 subtree through the coordinator, two arms at the same time on
# one machine (README.md). Usage, from anywhere in the repository:
#
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh base zone
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh warmy certx
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh zone2 base2   (replicates)
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh full ship2    (all levers; ship replicate)
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh fast wlowc    (Nov 6: the warm-call beta)
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh best bestcore (Nov 6: the core switches on top)
#   experiments/wil50_feasibility/pipeline_ab/run_pair.sh bestcore cpcut (Nov 7: the cut phase in CP form)
#   PAB_RUNS=.../runs_cpt experiments/wil50_feasibility/pipeline_ab/run_pair.sh bestcoret cpcutt
#                                        (Nov 7: the same in the production configuration, side 1)
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
#   fast   full + the Nov 5 speed switches (../speed/): R's rank term in fp32
#          off certificate iterations (DNN_R32_NONCERT=1) and the two-phase
#          beta for cold calls, 0.7 x beta for warm ones (DNN_BETA_SCHED)
#   cross  fast + the certificate at the predicted cutoff crossing between two
#          regular checks (DNN_ADMM_CROSS=2, ../speed/nov6_admm/): fathoms are
#          found 0-9 iterations earlier, nothing else changes
#   wlow   fast with the warm-call beta at PAB_WARM_W x base (default 0.25)
#          instead of 0.7, for children and cut rounds alike (../speed/nov6_admm/:
#          on 24 warm children -40% iterations to fathom, misses higher)
#   wlowc  wlow for the children, the cut rounds kept at 0.7 (the `wc` key):
#          the warm set measured children only
#   best   wlowc + DNN_ADMM_CROSS=2
#   core   fast + the three core switches of ../speed/nov6_core/ (the lazy
#          certificate DNN_CORE_LAZYCERT=1, no symmetrisation passes
#          DNN_CORE_NOSYMV=1 DNN_NOSYMW=1, the partial refresh
#          DNN_CORE_REFRESH_PARTIAL=1): the same decisions at about 1.1x less
#          CPU per iteration, so it is judged in DNN CPU, not iterations
#   bestcore  best + the same three core switches: every Nov 6 lever at once
#   wtwo   wlowc with a two-phase beta for a node's own warm solve
#          (DNN_ADMM_WARM_BETA=$PAB_WARM_TWO, default 0.25:60:0.7): 0.25 until the
#          first refresh at or after iteration 60 (80 at refresh 40), then 0.7, so
#          that fast children keep 0.25's gain and a miss hands its cut phase a
#          state like fast's; cut rounds stay at 0.7
#   cpcut  bestcore with the cut phase in CP form (QAP_DNN_CUT_CP=100:1000,
#          roadmap/admm/cpadmm/README.md "In the search"): one CP-ADMM solve with the
#          pool in the PRSM's Y-step and 1,000 new triangles separated every 100
#          iterations, the rounds' budget (12 x 200 + 1,000) and the stall rule per check
#   cpcut2 cpcut + a pooled child's first solve in CP form (QAP_DNN_POOL_CP=1)
#   cpcutl cpcut + the lazy certificate inside the CP solve (DNN_CP_LAZYCERT=1, with
#          bestcore's DNN_CORE_LAZYCERT=1): Rump's test only at checks that can decide,
#          the skipped ones certified at the exit with their own cut term
#   bestcoret, cpcutt  bestcore and cpcut in the production configuration of a real
#          wil50 run (../sym/README.md): the transposed instance (PAB_INSTANCE_T, default
#          qaplib/wil50t.dat) with its symmetry file (-s latest/SYM/wil50t.sym), the
#          regional cuts on side 1 (QAP_DNN_CUT_SIDE=1), and the node-10 spec transposed
#          (PAB_EXTRA_T, default -F 27:39,48:30,49:25 -x 0:4:1:300:1; i:j on wil50.dat is
#          j:i on wil50t.dat). There is no symmetry below node 10, so bestcoret's tree is
#          bestcore's up to the transposition's rounding. summarize.py flips their
#          traces' fixes back to wil50.dat's, so that their nodes match the other arms'
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

[ $# -ge 1 ] || { sed -n '2,/^set -u$/p' "$0" | sed '$d'; exit 2; }
for arm in "$@"; do
  case $arm in base|zone|warmy|certx|ship|full|zone2|base2|ship2|fast|cross|wlow|wlowc|best|core|bestcore|wtwo|cpcut|cpcut2|cpcutl|bestcoret|cpcutt) ;; *) echo "unknown arm '$arm' (base, zone, warmy, certx, ship, full, zone2, base2, ship2, fast, cross, wlow, wlowc, best, core, bestcore, wtwo, cpcut, cpcut2, cpcutl, bestcoret, cpcutt)"; exit 2;; esac
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
# the transposed arms (bestcoret, cpcutt): wil50t.dat, its symmetry file, side 1
INSTANCE_T=${PAB_INSTANCE_T:-qaplib/wil50t.dat}
EXTRA_T=${PAB_EXTRA_T:--F 27:39,48:30,49:25 -x 0:4:1:300:1 -s $R/latest/SYM/wil50t.sym}

# the arm whose settings a transposed arm takes
base_arm() { case $1 in bestcoret) echo bestcore;; cpcutt) echo cpcut;; *) echo "$1";; esac; }
is_transposed() { case $1 in bestcoret|cpcutt) return 0;; *) return 1;; esac; }

arm_env() {
  local arm=$1
  set -- "$(base_arm "$1")"
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
    fast)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7 ;;
    cross) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7 DNN_ADMM_CROSS=2 ;;
    wlow)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=${PAB_WARM_W:-0.25} ;;
    wlowc) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=${PAB_WARM_W:-0.25},wc=0.7 ;;
    best)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=${PAB_WARM_W:-0.25},wc=0.7 DNN_ADMM_CROSS=2 ;;
    core)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7 \
                  DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 ;;
    bestcore) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=${PAB_WARM_W:-0.25},wc=0.7 DNN_ADMM_CROSS=2 \
                  DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 ;;
    cpcut|cpcut2|cpcutl) export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=${PAB_WARM_W:-0.25},wc=0.7 DNN_ADMM_CROSS=2 \
                  DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 \
                  QAP_DNN_CUT_CP=${PAB_CUT_CP:-100:1000}
           if [ "$1" = cpcut2 ]; then export QAP_DNN_POOL_CP=1; fi
           if [ "$1" = cpcutl ]; then export DNN_CP_LAZYCERT=1; fi ;;
    wtwo)  export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 \
                  QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=${SHIP_BITS:-32} \
                  DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7 DNN_ADMM_WARM_BETA=${PAB_WARM_TWO:-0.25:60:0.7} ;;
    certx) export QAP_DNN_CUT_ZONE=0.005 DNN_CERTX=32 DNN_CERTX_USE=1 DNN_CERTX_QUIET=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=10 ;;
  esac
  if is_transposed "$arm"; then export QAP_DNN_CUT_SIDE=1; fi
}

run_arm() {
  local arm=$1 A=$RUNS/$1
  local B=local:$A/backend
  local base inst=$INSTANCE extra=$EXTRA
  base=$(base_arm "$arm")
  if is_transposed "$arm"; then inst=$INSTANCE_T; extra=$EXTRA_T; fi
  mkdir -p "$A"
  if [ -f "$A/backend/store/runs/$arm/results.txt" ]; then
    echo "$arm: finished earlier"
    return 0
  fi
  arm_env "$arm"
  if [ ! -x "$A/bin/bnbdrv" ]; then
    # a private copy, so that a rebuild during the run does not change it
    if { [ "$base" = warmy ] || [ "$base" = full ]; } && ! grep -aq QAP_DNN_CUT_DAMP1 latest/bnbdrv; then
      echo "warmy: this bnbdrv has no first-step damping yet: git pull, make solver, run again"; return 1
    fi
    if [ "$base" = fast ] && ! grep -aq DNN_BETA_SCHED latest/bnbdrv; then
      echo "fast: this bnbdrv has no DNN_BETA_SCHED or DNN_R32_NONCERT yet: git pull, make solver, run again"; return 1
    fi
    if { [ "$base" = cross ] || [ "$base" = best ] || [ "$base" = bestcore ]; } && ! grep -aq DNN_ADMM_CROSS latest/bnbdrv; then
      echo "$arm: this bnbdrv has no DNN_ADMM_CROSS yet: git pull, make solver, run again"; return 1
    fi
    if { [ "$base" = core ] || [ "$base" = bestcore ]; } && ! grep -aq DNN_CORE_LAZYCERT latest/bnbdrv; then
      echo "$arm: this bnbdrv has no DNN_CORE_* switches yet: git pull, make solver, run again"; return 1
    fi
    if { [ "$base" = wlowc ] || [ "$base" = best ] || [ "$base" = bestcore ]; } && ! grep -q "cut_round = 1" latest/dnncut.C; then
      # before the Nov 6 fix, wc also caught pooled children (they pass L_add), so wlowc was fast
      echo "$arm: this tree's wc key also applies to pooled children (fixed Nov 6): git pull, make solver, run again"; return 1
    fi
    if [ "$arm" = wtwo ] && ! grep -q "accel && !par.cut_round" latest/dnn.C; then
      echo "wtwo: this tree's DNN_ADMM_WARM_BETA also changes the cut rounds (fixed Nov 6): git pull, make solver, run again"; return 1
    fi
    if [ latest/bnbdrv -ot latest/dnn.C ] || [ latest/bnbdrv -ot latest/dnncut.C ] || [ latest/bnbdrv -ot latest/dnn.h ]; then
      echo "$arm: latest/bnbdrv is older than the DNN sources: make solver, run again"; return 1
    fi
    if { [ "$base" = cpcut ] || [ "$base" = cpcut2 ] || [ "$base" = cpcutl ]; } && ! grep -aq QAP_DNN_CUT_CP latest/bnbdrv; then
      echo "$arm: this bnbdrv has no CP cut phase (QAP_DNN_CUT_CP) yet: git pull, make solver, run again"; return 1
    fi
    if [ "$arm" = cpcutl ] && ! grep -aq DNN_CP_LAZYCERT latest/bnbdrv; then
      echo "cpcutl: this bnbdrv has no lazy certificate in CP form (DNN_CP_LAZYCERT) yet: git pull, make solver, run again"; return 1
    fi
    if [ "$arm" = cpcut2 ] && ! grep -aq QAP_DNN_POOL_CP latest/bnbdrv; then
      echo "cpcut2: this bnbdrv has no CP first solve for pooled children (QAP_DNN_POOL_CP) yet: git pull, make solver, run again"; return 1
    fi
    if { [ "$arm" = ship ] || [ "$arm" = ship2 ] || [ "$arm" = full ]; } && ! grep -aq QAP_DNN_SHIP latest/bnbdrv; then
      echo "ship: this bnbdrv cannot carry the DNN state yet: git pull, make solver, run again"; return 1
    fi
    if is_transposed "$arm" && { [ ! -f "$inst" ] || [ ! -f latest/SYM/wil50t.sym ]; }; then
      echo "$arm: $inst or latest/SYM/wil50t.sym is missing (../sym/make_sym.py writes both)"; return 1
    fi
    if [ "$arm" = cpcutt ] && ! grep -aq "regional cuts only (side 0 or 1)" latest/bnbdrv; then
      # before the side-1 projector (roadmap/admm/cpadmm, "Side 1") a side-1 CP phase silently ran the rounds
      echo "cpcutt: this bnbdrv's CP cut phase runs the rounds on side 1: git pull, make solver, run again"; return 1
    fi
    mkdir -p "$A/bin" && cp latest/bnbdrv "$A/bin/bnbdrv"
    echo "$(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)')" > "$A/bin/commit.txt"
  fi
  env | grep -E '^(QAP_|DNN_)' | sort > "$A/env.txt"
  echo "$(date -u +%FT%TZ) $CPU slots $SLOTS $(cat "$A/bin/commit.txt")" >> "$A/hosts.txt"
  if [ ! -f "$A/backend/ledger.db" ]; then
    $PY -m coord.qapctl init "$arm" --backend "$B" --instance "$inst" --params "$R/latest/PARAM/nug_deep_xbar9_all.par" \
      --ub "$UB" --slots "$SLOTS" --extra "$extra" --trace --solver "$A/bin/bnbdrv" \
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
