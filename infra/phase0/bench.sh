#!/bin/bash
# Phase 0 of docs/distributed/plan.md, the part that runs on the instance (README.md).
# Started by the user data that launch.sh writes, as root, from the extracted
# source tree, with QAP_BUCKET, QAP_OUT (the S3 prefix for results) and
# AWS_DEFAULT_REGION set. It builds the solver, runs the gates, the 32 wil50
# nodes of tuning/giveup/all32.txt one at a time on one core ("alone", 4 of
# them) and on every physical core at once ("loaded", all 32), uploads out/
# to S3 after every step, and shuts the instance down (which terminates it).
#
# Local test (no S3, no shutdown, any Linux): QAP_LOCAL=1 bash infra/phase0/bench.sh
#   QAP_ALONE="19"  the alone nodes (default "19 21 1 11")
#   QAP_LOADED="17 19 20 21"  the loaded nodes (default all 32)
#   QAP_CORES=2  at most this many cores for the loaded run
#   QAP_GATES="bruteforce"  the gates (default bruteforce resume parallel regression dnncheck)
set -uo pipefail
R=$(pwd)
OUT=$R/out
mkdir -p "$OUT"
exec > >(tee -a "$OUT/bench.log") 2>&1
LOCAL=${QAP_LOCAL:-0}
ALONE=${QAP_ALONE:-"19 21 1 11"}
LOADED=${QAP_LOADED:-$(seq 1 32 | tr '\n' ' ')}
GATES=${QAP_GATES:-"bruteforce resume parallel regression dnncheck"}
NODES=experiments/wil50_feasibility/tuning/giveup/all32.txt

up() {
  [ "$LOCAL" = 1 ] && return 0
  tar czf /tmp/out.tgz -C "$R" out && aws s3 cp --only-show-errors /tmp/out.tgz "s3://$QAP_BUCKET/$QAP_OUT/out.tgz"
}
step() { echo "$(date -u +%FT%TZ) $*" | tee -a "$OUT/steps.txt"; up; }
finish() {
  echo "$(date -u +%FT%TZ) end ($1)" >> "$OUT/steps.txt"
  up
  [ "$LOCAL" = 1 ] || shutdown -h now
}
trap 'finish exit' EXIT

step "start"
# --- the machine
if [ "$LOCAL" != 1 ]; then
  T=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H 'X-aws-ec2-metadata-token-ttl-seconds: 3600')
  md() { curl -s -H "X-aws-ec2-metadata-token: $T" "http://169.254.169.254/latest/meta-data/$1"; }
  { echo "instance_type $(md instance-type)"; echo "az $(md placement/availability-zone)"
    echo "lifecycle $(md instance-life-cycle)"; echo "instance_id $(md instance-id)"; } > "$OUT/sysinfo.txt"
else
  echo "instance_type local-$(hostname)" > "$OUT/sysinfo.txt"
fi
CORES=$(lscpu -p=CPU,CORE,SOCKET | grep -v '^#' | sort -t, -k3,3n -k2,2n -k1,1n | awk -F, '!seen[$3","$2]++ {print $1}')
[ -n "${QAP_CORES:-}" ] && CORES=$(echo "$CORES" | head -n "$QAP_CORES")
NC=$(echo "$CORES" | wc -l)
{ echo "physical_cores $NC"; echo "vcpus $(nproc)"; echo "arch $(uname -m)"
  echo "model $(lscpu | sed -n 's/^Model name: *//p' | head -1)"
  echo "amx_bf16 $(grep -m1 -c amx_bf16 /proc/cpuinfo)"
  echo "mem_gb $(free -g | awk '/^Mem/ {print $2}')"; echo "kernel $(uname -r)"; } >> "$OUT/sysinfo.txt"
lscpu > "$OUT/lscpu.txt"
step "machine: $(tr '\n' ' ' < "$OUT/sysinfo.txt")"

# --- the toolchain (Amazon Linux 2023; the local test uses what is there)
PY=python3
if [ "$LOCAL" != 1 ]; then
  dnf install -y -q gcc-c++ make python3.11 python3.11-pip tar gzip util-linux > "$OUT/dnf.log" 2>&1 || { step "dnf failed"; exit 1; }
  PY=python3.11
  $PY -m pip install -q numpy scipy-openblas32 > "$OUT/pip.log" 2>&1 || { step "pip failed"; exit 1; }
fi
{ g++ --version | head -1; echo "cc: $(cc --version | head -1)"; $PY --version
  $PY -c "import scipy_openblas32 as o; print('scipy-openblas32', o.get_lib_dir())"; } > "$OUT/toolchain.txt" 2>&1
step "toolchain: $(tr '\n' ' ' < "$OUT/toolchain.txt")"

# --- the build
if ! { make solver PYTHON=$PY && make dnntest PYTHON=$PY; } > "$OUT/build.log" 2>&1; then
  step "build FAILED (out/build.log)"
  exit 1
fi
ldd latest/dnntest | grep -iE 'blas|lapack' > "$OUT/blas.txt"
step "build ok: $(tr '\n' ' ' < "$OUT/blas.txt")"

# --- the gates (results are recorded; a failure does not stop the benchmark)
: > "$OUT/gates.txt"
for g in $GATES; do
  t0=$(date +%s)
  if make "$g" PYTHON=$PY > "$OUT/gate_$g.log" 2>&1; then r=pass; else r=FAIL; fi
  echo "$g $r $(( $(date +%s) - t0 ))s" >> "$OUT/gates.txt"
done
step "gates: $(tr '\n' ' ' < "$OUT/gates.txt")"

# --- the benchmark: production DNN flags (bundle D), AMX where the CPU has it
run_node() {   # node line, cpu, directory, [extra flags]
  local k=$1 cpu=$2 dir=$3 spec
  spec=$(sed -n "${k}p" "$NODES" | awk '{print $1}')
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=1 \
    taskset -c "$cpu" latest/dnntest -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i "${4:-2500}" -u 48816 -v \
    -x "$spec" qaplib/wil50.dat > "$dir/node_$k.log" 2>&1
}
first=$(echo "$CORES" | head -1)
mkdir -p "$OUT/alone" "$OUT/loaded"
for k in $ALONE; do
  run_node "$k" "$first" "$OUT/alone"
  step "alone node $k: $(grep -o 'status=[a-z]* iters=[0-9]* cpu=[0-9.]* ms_per_iter=[0-9.]*' "$OUT/alone/node_$k.log")"
done

# loaded: one process per physical core, a shared queue; once it is empty the
# idle cores run short filler solves, so the load stays full until the last
# real node ends
Q=$OUT/loaded/queue
printf '%s\n' $LOADED > "$Q"
: > "$OUT/loaded/finished"
rm -f "$OUT/loaded/done"
worker() {
  local cpu=$1 k
  mkdir -p "$OUT/loaded/filler_$cpu"
  while :; do
    k=$(flock "$Q.lock" sh -c "head -1 '$Q'; sed -i 1d '$Q'")
    [ -z "$k" ] && break
    run_node "$k" "$cpu" "$OUT/loaded"
    echo "$k" >> "$OUT/loaded/finished"
  done
  while [ ! -f "$OUT/loaded/done" ]; do run_node 1 "$cpu" "$OUT/loaded/filler_$cpu" 40; done
  rm -rf "$OUT/loaded/filler_$cpu"
}
n_loaded=$(echo $LOADED | wc -w)
t0=$(date +%s)
pids=()
for c in $CORES; do worker "$c" & pids+=($!); done
while [ "$(wc -l < "$OUT/loaded/finished")" -lt "$n_loaded" ]; do sleep 10; done
echo "loaded_wall_s $(( $(date +%s) - t0 ))" >> "$OUT/sysinfo.txt"
touch "$OUT/loaded/done"
for p in "${pids[@]}"; do wait "$p"; done
step "loaded: $n_loaded nodes on $NC cores in $(( $(date +%s) - t0 )) s"

$PY infra/phase0/summarize.py instance "$OUT" > "$OUT/summary.txt" 2>&1
step "summary: $(head -c 400 "$OUT/summary.txt" | tr '\n' ' ')"
