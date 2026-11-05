#!/bin/bash
# Sequential A/B on wil50 nodes of tuning/giveup/all32.txt, pinned to core 3,
# one process at a time, production-like flags (bundle D without AMX):
# PRSM -g 1.0 -G 1.0, degree 3, -L 2, -r 40, -H 300, -i 2500 -u 48816,
# DNN_CHOLQR_PASSES=1, one BLAS thread. Copied from
# ../../solvers/accel/run_all32.sh with -P 1 and taskset -c 3.
# usage: run_nodes.sh name "ENV=1 ..." "id id ..." [dnntest]
# Logs: runs/<name>_<id>.log (id = line number in all32.txt); resumable.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
name=$1; envs=$2; ids=$3; bin=${4:-$R/latest/dnntest}
mkdir -p "$D/runs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
for k in $ids; do
  spec=$(sed -n "${k}p" "$D/../../tuning/giveup/all32.txt")
  f=$D/runs/${name}_$k.log
  grep -q '^RESULT' "$f" 2>/dev/null && continue
  env $envs taskset -c 3 "$bin" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 2500 -u 48816 -v -v \
    -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
done
