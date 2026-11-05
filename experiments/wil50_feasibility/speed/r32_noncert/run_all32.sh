#!/bin/bash
# A/B of DNN_R32_NONCERT=1 (R's rank term in fp32 except on certificate
# iterations) on the 32 wil50 nodes of tuning/giveup/all32.txt, adapted from
# ../../solvers/accel/run_all32.sh: production-like flags (degree 3, -L 2,
# PRSM -g 1.0 -G 1.0, -r 40, -H 300, DNN_CHOLQR_PASSES=1; no AMX), one
# process at a time pinned to core 0 (-P 1), single-threaded BLAS. The
# configurations run interleaved per node (base, then nc, on the same node),
# so a change in the host's background load hits both alike.
#   base  the default R-step (fp64 dsyrk every iteration)
#   nc    DNN_R32_NONCERT=1
# NODES: line numbers of all32.txt to run (default all 32), e.g. NODES="1 5 9".
# CORE: the core to pin to (default 0). Resumable: a node whose log has a
# RESULT line is skipped. Logs: runs/<name>_<id>.log, id = line in all32.txt.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CONFIGS=${CONFIGS:-"base:-d_3_-L_2_-r_40: nc:-d_3_-L_2_-r_40:DNN_R32_NONCERT=1"}
NODES=${NODES:-$(seq 1 32)}
CORE=${CORE:-0}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
one() {
  name=$1; flags=${2//_/ }; envs=$3; k=$4; spec=$5
  f=$D/runs/${name}_$k.log
  grep -q '^RESULT' "$f" 2>/dev/null && return 0
  env $envs taskset -c "$CORE" "$R/latest/dnntest" -T 1 $flags -g 1.0 -G 1.0 -H 300 -i 2500 \
    -u 48816 -v -v -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
mkdir -p "$D/runs"
for k in $NODES; do
  spec=$(sed -n "${k}p" "$D/../../tuning/giveup/all32.txt" | awk '{print $1}')
  for c in $CONFIGS; do
    IFS=: read -r name flags envs <<< "$c"
    one "$name" "$flags" "$envs" "$k" "$spec"
  done
done
