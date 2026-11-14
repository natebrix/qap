#!/bin/bash
# Phase 1 of docs/distributed/plan.md, the part that runs on the instance
# (README.md). Started by the user data that launch.sh writes, as root, from
# the extracted source tree.
#
# QAP_ROLE=all (the default): build the solver, run the gates, create the run
# (`qapctl init`: the frontier pass, here, so the whole tree comes from one
# build) unless it exists, then run the coordinator and one worker with a slot
# per physical core until the run is finished. Then, with QAP_REFERENCE=1,
# bench/parallel.py with the same settings on the same cores, and compare.py:
# the acceptance check. A run that exists already (a relaunch after a Spot
# interruption) is resumed: the coordinator replays what was merged, the
# worker takes the queued units.
# QAP_ROLE=worker: the build, the gates and a worker, for a run that another
# instance coordinates.
#
# A Spot interruption notice drains the worker: every unit goes back with its
# open nodes (or its last checkpoint) before the instance is reclaimed.
# out/ goes to the store as runs/RUN/hosts/HOST/out.tgz after every step and
# every 5 minutes. At the end the instance shuts down, which terminates it.
#
# Set by the user data: QAP_RUN, QAP_BACKEND (aws:BUCKET?region=R), QAP_BUCKET
# (for uploads before boto3 is installed), QAP_RECIPE (a name in recipes/;
# needed only to create the run), QAP_ROLE, QAP_NO_UB=1 (no incumbent),
# QAP_REFERENCE (default 1), QAP_GATES (default below; "none" skips them),
# QAP_SLOTS (default: the physical cores), QAP_VISIBILITY (public or private:
# the coordinator's public status, coord/publish.py; for a new run only).
#
# Local test (any Linux; no dnf, no shutdown): QAP_LOCAL=1, QAP_BACKEND pointing
# at a moto server (aws:B?region=eu-north-1&endpoint=http://127.0.0.1:PORT),
# QAP_HOST to tell two local nodes apart, and QAP_IMDS=http://127.0.0.1:PORT2 for
# a stand-in of the instance metadata service (without it the Spot watcher does
# not run in local mode).
set -uo pipefail
R=$(pwd)
OUT=$R/out
mkdir -p "$OUT" "$R/tmp" "$R/work"
exec > >(tee -a "$OUT/node.log") 2>&1
export TMPDIR=$R/tmp HOME=${HOME:-/root}   # /tmp is a small tmpfs on Amazon Linux 2023
LOCAL=${QAP_LOCAL:-0}
RUN=${QAP_RUN:?QAP_RUN is not set}
BK=${QAP_BACKEND:?QAP_BACKEND is not set}
ROLE=${QAP_ROLE:-all}
NO_UB=${QAP_NO_UB:-0}
REFERENCE=${QAP_REFERENCE:-1}
GATES=${QAP_GATES:-"bruteforce resume parallel regression dnncheck coord"}
[ "$GATES" = none ] && GATES=""
IMDS=${QAP_IMDS:-}
[ "$LOCAL" = 1 ] || IMDS=${IMDS:-http://169.254.169.254}
PY=python3
HAVE_BK=0
case $ROLE in all|worker) ;; *) echo "QAP_ROLE must be all or worker"; exit 2 ;; esac

imds() {   # GET a metadata path (IMDSv2); empty when there is no metadata service
  local tok
  [ -n "$IMDS" ] || return 0
  tok=$(curl -s -m 2 -X PUT "$IMDS/latest/api/token" -H 'X-aws-ec2-metadata-token-ttl-seconds: 60') || return 0
  curl -s -m 2 -f -H "X-aws-ec2-metadata-token: $tok" "$IMDS/latest/meta-data/$1" || true
}
HOST=${QAP_HOST:-}
[ -n "$HOST" ] || HOST=$(imds instance-id)
[ -n "$HOST" ] || HOST=$(hostname)
bk() { $PY infra/phase1/bk.py "$BK" "$@"; }
up() {   # out/ to the store; never fatal
  tar czf "$TMPDIR/out.tgz" -C "$R" out 2> /dev/null || return 0
  if [ "$HAVE_BK" = 1 ]; then bk put "runs/$RUN/hosts/$HOST/out.tgz" "$TMPDIR/out.tgz" || echo "upload failed"
  elif [ "$LOCAL" != 1 ] && [ -n "${QAP_BUCKET:-}" ]; then
    aws s3 cp --only-show-errors "$TMPDIR/out.tgz" "s3://$QAP_BUCKET/runs/$RUN/hosts/$HOST/out.tgz" || echo "upload failed"
  fi
}
step() { echo "$(date -u +%FT%TZ) $*" | tee -a "$OUT/steps.txt"; up; }
finish() {
  echo "$(date -u +%FT%TZ) end ($1)" >> "$OUT/steps.txt"
  touch "$OUT/finished"
  for f in worker coordinator; do [ -f "$OUT/$f.pid" ] && kill -TERM "$(cat "$OUT/$f.pid")" 2> /dev/null; done
  kill $(jobs -p) 2> /dev/null
  up
  [ "$LOCAL" = 1 ] || shutdown -h now
}
trap 'finish exit' EXIT

step "start: run $RUN, role $ROLE, host $HOST"
# --- the machine
{ echo "instance_type $(imds instance-type)"; echo "az $(imds placement/availability-zone)"
  echo "lifecycle $(imds instance-life-cycle)"; echo "host $HOST"; } > "$OUT/sysinfo.txt"
CORES=$(lscpu -p=CPU,CORE,SOCKET | grep -v '^#' | sort -t, -k3,3n -k2,2n -k1,1n | awk -F, '!seen[$3","$2]++ {print $1}')
NC=$(echo "$CORES" | wc -l)
SLOTS=${QAP_SLOTS:-$NC}
{ echo "physical_cores $NC"; echo "slots $SLOTS"; echo "vcpus $(nproc)"; echo "arch $(uname -m)"
  echo "model $(lscpu | sed -n 's/^Model name: *//p' | head -1)"
  echo "amx_bf16 $(grep -m1 -c amx_bf16 /proc/cpuinfo)"
  echo "mem_gb $(free -g | awk '/^Mem/ {print $2}')"; echo "kernel $(uname -r)"; } >> "$OUT/sysinfo.txt"
step "machine: $(tr '\n' ' ' < "$OUT/sysinfo.txt")"

# OpenBLAS on Graviton3 (Neoverse V1): its auto-selected kernels give the DNN certificate
# -inf at m >= 9 (Nov 13, scipy-openblas32 on c7g.xlarge: nug12's root -inf against the
# prototype's 567.98; tests/dnn_test.py fails, the parallel gate's DNN check never ends).
# The Neoverse N1 kernels give x86's values to the digit. --node-env OPENBLAS_CORETYPE=...
# overrides this; OPENBLAS_CORETYPE=auto leaves the choice to OpenBLAS (a diagnosis).
if [ "${OPENBLAS_CORETYPE:-}" = auto ]; then unset OPENBLAS_CORETYPE
elif [ "$(uname -m)" = aarch64 ] && [ -z "${OPENBLAS_CORETYPE:-}" ]; then export OPENBLAS_CORETYPE=NEOVERSEN1; fi
echo "openblas_coretype ${OPENBLAS_CORETYPE:-auto}" >> "$OUT/sysinfo.txt"

# --- the toolchain (Amazon Linux 2023, as in infra/phase0/bench.sh, plus boto3)
if [ "$LOCAL" != 1 ]; then
  dnf install -y -q gcc-c++ make python3.11 python3.11-pip tar gzip util-linux > "$OUT/dnf.log" 2>&1 || { step "dnf failed"; exit 1; }
  PY=python3.11
  $PY -m pip install -q numpy scipy-openblas32 boto3 > "$OUT/pip.log" 2>&1 || { step "pip failed"; exit 1; }
fi
PY=${QAP_PYTHON:-$PY}
$PY -c "import boto3" 2> /dev/null && HAVE_BK=1
{ g++ --version | head -1; echo "cc: $(cc --version | head -1)"; $PY --version
  $PY -c "import scipy_openblas32 as o; print('scipy-openblas32', o.get_lib_dir())" 2> /dev/null || echo "no scipy-openblas32 (system LAPACK)"
  $PY -c "import boto3; print('boto3', boto3.__version__)" 2> /dev/null || echo "no boto3"; } > "$OUT/toolchain.txt" 2>&1
step "toolchain: $(tr '\n' ' ' < "$OUT/toolchain.txt")"
[ "$HAVE_BK" = 1 ] || { step "no boto3 for $PY"; exit 1; }
# from here on out/ also goes up every 5 minutes, whatever the main line is doing
( while [ ! -f "$OUT/finished" ]; do sleep 300; up > /dev/null; done ) &

# --- the build and the gates (a failed gate stops the node: no units on a build that fails them)
if ! make solver PYTHON=$PY > "$OUT/build.log" 2>&1; then step "build FAILED (out/build.log)"; exit 1; fi
ldd latest/bnbdrv | grep -iE 'blas|lapack' > "$OUT/blas.txt"
step "build ok: ${QAP_COMMIT:-$(git -C "$R" rev-parse --short HEAD 2> /dev/null)} $(tr '\n' ' ' < "$OUT/blas.txt")"
: > "$OUT/gates.txt"
# QAP_PROBE (--node-env on launch.sh): a shell command run in latest/ after the build,
# its output in out/probe.log, before the gates (a diagnosis on the instance's CPU
# without logging in; Nov 13: the DNN certificate on Graviton)
if [ -n "${QAP_PROBE:-}" ]; then
  ( cd latest && timeout 1200 bash -c "$QAP_PROBE" ) > "$OUT/probe.log" 2>&1; prc=$?
  step "probe exit $prc: $(wc -l < "$OUT/probe.log") lines in out/probe.log"
  up > /dev/null
fi
# python -u: a gate's progress reaches its log (and the 5-minute upload) while it runs;
# timeout: a hung gate is a failure after 30 minutes, not an idle instance until --limit
# (Nov 13: tests/parallel_test.py hung on the first real c7g.xlarge with an empty log)
for g in $GATES; do
  t0=$(date +%s)
  timeout 1800 make "$g" PYTHON="$PY -u" > "$OUT/gate_$g.log" 2>&1; rc=$?
  case $rc in 0) r=pass ;; 124) r=TIMEOUT ;; *) r=FAIL ;; esac
  echo "$g $r $(( $(date +%s) - t0 ))s" >> "$OUT/gates.txt"
done
step "gates: $(tr '\n' ' ' < "$OUT/gates.txt")"
grep -q -E 'FAIL|TIMEOUT' "$OUT/gates.txt" && { step "a gate failed: stopping (out/gate_*.log)"; exit 1; }

# --- the Spot watcher: an interruption notice drains the worker (2 minutes' warning)
spot_watch() {
  local tok code
  [ -n "$IMDS" ] || return 0
  while [ ! -f "$OUT/finished" ]; do
    sleep "${QAP_SPOT_POLL:-5}"
    tok=$(curl -s -m 2 -X PUT "$IMDS/latest/api/token" -H 'X-aws-ec2-metadata-token-ttl-seconds: 60') || continue
    code=$(curl -s -m 2 -o "$OUT/spot_action.json" -w '%{http_code}' -H "X-aws-ec2-metadata-token: $tok" \
      "$IMDS/latest/meta-data/spot/instance-action") || continue
    [ "$code" = 200 ] || continue
    touch "$OUT/draining"
    step "Spot interruption notice $(cat "$OUT/spot_action.json"): draining the worker"
    [ -f "$OUT/worker.pid" ] && kill -TERM "$(cat "$OUT/worker.pid")" 2> /dev/null
    return 0
  done
}
spot_watch &

# --- the worker and the coordinator, each restarted if it dies before the run is finished
worker_loop() {
  local n=0 rc
  while [ ! -f "$OUT/draining" ] && [ ! -f "$OUT/stopping" ]; do
    $PY -m coord.qapctl worker "$RUN" --backend "$BK" --slots "$SLOTS" --workdir "$R/work" >> "$OUT/worker.log" 2>&1 &
    echo $! > "$OUT/worker.pid"
    wait $!; rc=$?
    rm -f "$OUT/worker.pid"
    [ "$(bk status "$RUN")" = finished ] && return 0
    [ -f "$OUT/draining" ] || [ -f "$OUT/stopping" ] && return 1
    n=$((n + 1)); step "worker exited ($rc) before the run finished; restart $n of 5"
    [ $n -ge 5 ] && return 1
    sleep 10
  done
  return 1
}
coord_loop() {
  local n=0 rc
  while :; do
    $PY -m coord.qapctl coordinator "$RUN" --backend "$BK" >> "$OUT/coordinator.log" 2>&1 &
    echo $! > "$OUT/coordinator.pid"
    wait $!; rc=$?
    rm -f "$OUT/coordinator.pid"
    [ $rc = 0 ] && return 0
    [ -f "$OUT/draining" ] && return 1
    n=$((n + 1)); step "coordinator exited ($rc); restart $n of 5 (it replays what was merged)"
    if [ $n -ge 5 ]; then   # nothing would merge: stop the worker too (its units go back)
      touch "$OUT/stopping"
      [ -f "$OUT/worker.pid" ] && kill -TERM "$(cat "$OUT/worker.pid")" 2> /dev/null
      # and the instance: the marker above needs a disk with room, and tho40-a's
      # was full (Nov 14), so its worker loop never saw it and the instance idled
      # for hours; finish (out/ uploaded, then shutdown) needs none
      [ "$LOCAL" = 1 ] || finish "the coordinator failed 5 times"
      return 1
    fi
    sleep 10
  done
}

ST=$(bk status "$RUN") || { step "cannot read run $RUN from $BK: stopping"; exit 1; }
if [ "$ROLE" = worker ]; then
  for i in $(seq 60); do [ "$ST" != none ] && break; sleep 30; ST=$(bk status "$RUN"); done
  [ "$ST" = none ] && { step "no run $RUN after 30 minutes: stopping"; exit 1; }
  [ "$ST" = finished ] && { step "run $RUN is finished: nothing to do"; exit 0; }
  step "worker: $SLOTS slots on run $RUN ($ST)"
  worker_loop; rc=$?
  step "worker: $([ $rc = 0 ] && echo "run $RUN finished" || echo stopped) ($(grep -c ' -> ' "$OUT/worker.log" 2> /dev/null) units)"
  exit 0
fi

# role all: the run
if [ "$ST" = finished ]; then step "run $RUN is finished already: nothing to do"; exit 0; fi
if [ "$ST" = none ]; then
  RF=infra/phase1/recipes/${QAP_RECIPE:-}.sh
  [ -n "${QAP_RECIPE:-}" ] && [ -f "$RF" ] || { step "no run $RUN and no recipe (QAP_RECIPE=${QAP_RECIPE:-}): stopping"; exit 1; }
  INSTANCE= PARAMS= OPT= UB= EXTRA= ENV= COORD=
  source "$RF"
  [ "$NO_UB" = 1 ] && UB=
  # the recipe as this run uses it, for the reference run and for a relaunch
  { echo "# run $RUN: recipes/$QAP_RECIPE.sh$([ "$NO_UB" = 1 ] && echo ', --no-ub')"; declare -p INSTANCE PARAMS OPT UB EXTRA ENV COORD; } > "$OUT/recipe.sh"
  args=(--instance "$INSTANCE" --params "$PARAMS" --slots "$SLOTS" "--extra=$EXTRA")
  for kv in $ENV; do args+=(--env "$kv"); done
  [ -n "$UB" ] && args+=(--ub "$UB")
  args+=(--visibility "${QAP_VISIBILITY:-public}" $COORD)
  step "init: run $RUN from recipes/$QAP_RECIPE.sh$([ -n "$UB" ] && echo ", incumbent $UB" || echo ", no incumbent"): the frontier pass"
  if ! $PY -m coord.qapctl init "$RUN" --backend "$BK" "${args[@]}" > "$OUT/init.log" 2>&1; then
    step "init FAILED (out/init.log)"; exit 1
  fi
  bk put "runs/$RUN/recipe.sh" "$OUT/recipe.sh"
  step "init: $(grep -o 'frontier at depth .*' "$OUT/init.log")"
else
  bk get "runs/$RUN/recipe.sh" "$OUT/recipe.sh" || : > "$OUT/recipe.sh"
  step "run $RUN exists ($ST): resuming"
fi
coord_loop & CL=$!
step "coordinator and worker ($SLOTS slots) running"
worker_loop; WR=$?
if [ -f "$OUT/draining" ]; then
  [ -f "$OUT/coordinator.pid" ] && kill -TERM "$(cat "$OUT/coordinator.pid")" 2> /dev/null
  wait $CL
  step "drained ($(grep -c ' -> ' "$OUT/worker.log" 2> /dev/null) units here): relaunch with infra/phase1/launch.sh $RUN"
  exit 0
fi
wait $CL; CR=$?
[ -f "$OUT/stopping" ] && { step "stopped: the coordinator failed 5 times (out/coordinator.log)"; exit 1; }
if [ $WR != 0 ] || [ $CR != 0 ] || ! bk get "runs/$RUN/results.txt" "$OUT/results.txt"; then
  step "the run did not finish (worker $WR, coordinator $CR): see out/worker.log, out/coordinator.log"; exit 1
fi
step "finished: $(grep -E 'Nodes:|best objective' "$OUT/results.txt" | tr -s ' ' | tr '\n' ' ')"

# --- the acceptance check: bench/parallel.py with the same settings, here
INSTANCE= PARAMS= OPT= UB= EXTRA= ENV= COORD=
source "$OUT/recipe.sh"
cmp=(); [ -n "$OPT" ] && cmp+=(--opt "$OPT")
if [ "$REFERENCE" = 1 ] && [ -n "$INSTANCE" ]; then
  ref=(); [ -n "$UB" ] && ref+=(--ub "$UB") && cmp+=(--fixed-ub)
  step "reference: bench/parallel.py --jobs $SLOTS on the same settings"
  env $ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PY bench/parallel.py "$INSTANCE" --jobs "$SLOTS" "${ref[@]}" \
    --params "$PARAMS" "--extra=$EXTRA" > "$OUT/reference.log" 2>&1 || step "reference run failed (out/reference.log)"
  $PY infra/phase1/compare.py "$OUT/results.txt" "$OUT/reference.log" "${cmp[@]}" > "$OUT/compare.txt"
else
  $PY infra/phase1/compare.py "$OUT/results.txt" "${cmp[@]}" > "$OUT/compare.txt"
fi
bk put "runs/$RUN/compare.txt" "$OUT/compare.txt"
step "check: $(sed -n '/^PASS\|^FAIL/,$p' "$OUT/compare.txt" | tr -s ' \n' ' ')"
