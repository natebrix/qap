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
#   QAP_INSTALL=1  also install the toolchain as on EC2 (Amazon Linux 2023, as root): a host
#                  that launch.sh did not start, e.g. a Lightsail instance (README.md)
#   QAP_ALONE="19"  the alone nodes (default "19 21 1 11")
#   QAP_LOADED="17 19 20 21"  the loaded nodes (default all 32)
#   QAP_CORES=2  at most this many cores for the loaded run
#   QAP_GATES="bruteforce"  the gates (default bruteforce resume parallel regression dnncheck)
#   QAP_BLAS="openblas openblas-skylakex blis"  BLAS variants, each benchmarked in out/blas_NAME/ (README.md)
set -uo pipefail
R=$(pwd)
OUT=$R/out
mkdir -p "$OUT"
exec > >(tee -a "$OUT/bench.log") 2>&1
LOCAL=${QAP_LOCAL:-0}
INSTALL=${QAP_INSTALL:-$((1 - LOCAL))}
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

# --- the toolchain (Amazon Linux 2023; the local test uses what is there unless QAP_INSTALL=1)
PY=python3
if [ "$INSTALL" = 1 ]; then
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

# --- BLAS variants (QAP_BLAS, README.md "BLAS variants"): each prepares its own dnntest and environment.
# Unset: the build above, as before. A variant that cannot be prepared is recorded and skipped.
DNNTEST=latest/dnntest; BENV=()
corename() {   # the kernels the wheel's OpenBLAS picks under the current BENV
  env "${BENV[@]}" $PY -c "import ctypes, glob, scipy_openblas32 as o
l = ctypes.CDLL(glob.glob(o.get_lib_dir() + '/libscipy_openblas*.so')[0]); f = l.scipy_openblas_get_corename; f.restype = ctypes.c_char_p; print(f().decode())" 2> /dev/null || echo "?"
}
build_variant() {   # name, make variables...: a copy of latest/ built on its own (dnn.o depends on the BLAS defines)
  local v=$1; shift
  rm -rf "latest_$v" && cp -r latest "latest_$v" && make -C "latest_$v" clean > /dev/null 2>&1
  make -C "latest_$v" dnntest PYTHON=$PY "$@" > "$OUT/build_$v.log" 2>&1 || return 1
  DNNTEST=latest_$v/dnntest
}
gcc14() {   # Amazon Linux 2023's gcc 14 (it knows znver4; the default gcc 11 does not)
  [ "$INSTALL" = 1 ] && dnf install -y -q gcc14 gcc14-c++ >> "$OUT/dnf.log" 2>&1
  G14=$(command -v gcc14-g++ || command -v g++-14) && C14=$(command -v gcc14-gcc || command -v gcc-14)
}
BLIS_TAG=${QAP_BLIS_TAG:-5.1}   # github.com/amd/blis (AOCL-BLAS), built here for this CPU
blis() {   # /opt/blis: one process per core is the solver's parallelism, so no threading
  [ -f /opt/blis/lib/libblis.so ] && return 0
  [ "$INSTALL" = 1 ] || return 1
  dnf install -y -q git lapack >> "$OUT/dnf.log" 2>&1 || return 1
  gcc14 || C14=gcc
  { rm -rf /opt/blis-src && git clone -q --depth 1 -b "$BLIS_TAG" https://github.com/amd/blis /opt/blis-src &&
    cd /opt/blis-src && CC=$C14 ./configure --prefix=/opt/blis --disable-threading --disable-static auto &&
    make -j "$(nproc)" && make install; } > "$OUT/build_blis_lib.log" 2>&1
  local rc=$?; cd "$R"; return $rc
}
prep_variant() {
  local v=$1
  DNNTEST=latest/dnntest; BENV=()
  case $v in
    openblas) ;;                                                    # the wheel, its own choice of kernels
    openblas-*) BENV=("OPENBLAS_CORETYPE=$(echo "${v#openblas-}" | tr a-z A-Z)") ;;   # the wheel, these kernels
    gcc14) gcc14 && build_variant "$v" CPP="$G14" CC="$C14" ;;   # the wheel; the solver compiled by gcc 14
    netlib) build_variant "$v" DNN_LAPACK="-l:liblapack.so.3 -l:libblas.so.3" DNN_LAPACK_DEFS= ;;   # reference BLAS: a check of the mechanism
    blis) blis && build_variant "$v" DNN_LAPACK="-L/opt/blis/lib -lblis -l:liblapack.so.3 -Wl,-rpath,/opt/blis/lib" DNN_LAPACK_DEFS= ;;
    *) return 1 ;;
  esac
}

# --- the benchmark: production DNN flags (bundle D), AMX where the CPU has it
run_node() {   # node line, cpu, directory, [extra flags]
  local k=$1 cpu=$2 dir=$3 spec
  spec=$(sed -n "${k}p" "$NODES" | awk '{print $1}')
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "${BENV[@]}" \
    taskset -c "$cpu" "$DNNTEST" -p infra/phase0/dnntest.par -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i "${4:-2500}" -u 48816 -v \
    -x "$spec" qaplib/wil50.dat > "$dir/node_$k.log" 2>&1
}
first=$(echo "$CORES" | head -1)
n_loaded=$(echo $LOADED | wc -w)
worker() {   # cpu, directory
  local cpu=$1 O=$2 k
  mkdir -p "$O/loaded/filler_$cpu"
  while :; do
    k=$(flock "$O/loaded/queue.lock" sh -c "head -1 '$O/loaded/queue'; sed -i 1d '$O/loaded/queue'")
    [ -z "$k" ] && break
    run_node "$k" "$cpu" "$O/loaded"
    echo "$k" >> "$O/loaded/finished"
  done
  while [ ! -f "$O/loaded/done" ]; do run_node 1 "$cpu" "$O/loaded/filler_$cpu" 40; done
  rm -rf "$O/loaded/filler_$cpu"
}
bench() {   # directory, label: alone, then loaded, then the summary
  local O=$1 lab=$2 k t0 c p pids=()
  mkdir -p "$O/alone" "$O/loaded"
  for k in $ALONE; do
    run_node "$k" "$first" "$O/alone"
    step "${lab}alone node $k: $(grep -o 'status=[a-z]* iters=[0-9]* cpu=[0-9.]* ms_per_iter=[0-9.]*' "$O/alone/node_$k.log")"
  done
  # loaded: one process per physical core, a shared queue; once it is empty the
  # idle cores run short filler solves, so the load stays full until the last
  # real node ends
  printf '%s\n' $LOADED > "$O/loaded/queue"
  : > "$O/loaded/finished"
  rm -f "$O/loaded/done"
  t0=$(date +%s)
  for c in $CORES; do worker "$c" "$O" & pids+=($!); done
  while [ "$(wc -l < "$O/loaded/finished")" -lt "$n_loaded" ]; do sleep 10; done
  echo "loaded_wall_s $(( $(date +%s) - t0 ))" >> "$O/sysinfo.txt"
  touch "$O/loaded/done"
  for p in "${pids[@]}"; do wait "$p"; done
  step "${lab}loaded: $n_loaded nodes on $NC cores in $(( $(date +%s) - t0 )) s"
  $PY infra/phase0/summarize.py instance "$O" > "$O/summary.txt" 2>&1
  step "${lab}summary: $(head -c 400 "$O/summary.txt" | tr '\n' ' ')"
}
if [ -z "${QAP_BLAS:-}" ]; then
  bench "$OUT" ""
else
  for v in $QAP_BLAS; do
    if ! prep_variant "$v"; then step "blas $v: not prepared (out/build_$v.log, out/build_blis_lib.log, out/dnf.log)"; continue; fi
    O=$OUT/blas_$v; mkdir -p "$O"
    cp "$OUT/sysinfo.txt" "$OUT/gates.txt" "$O/"
    { echo "blas $v"; echo "blas_core $(ldd "$DNNTEST" | grep -q scipy_openblas && corename || echo -)"; echo "blas_libs $(ldd "$DNNTEST" | grep -iE 'blas|lapack|blis' | awk '{print $1}' | tr '\n' ' ')"; } >> "$O/sysinfo.txt"
    step "blas $v: $DNNTEST ${BENV[*]:-} | $(tail -2 "$O/sysinfo.txt" | tr '\n' ' ')"
    bench "$O" "$v "
  done
fi
