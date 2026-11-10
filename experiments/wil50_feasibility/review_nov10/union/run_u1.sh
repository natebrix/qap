#!/bin/bash
# usage: run_u1.sh CPU NAME PARENT [dnntest args...]   (log: u1/NAME.log)
R=$(cd "$(dirname "$0")/../../../.." && pwd)
U=$(cd "$(dirname "$0")" && pwd)   # logs here; dnntest_u1 = latest/dnntest built with union_pass.patch
CPU=$1; NAME=$2; PARENT=$3; shift 3
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c "$CPU" "$U/dnntest_u1" \
  -p "$R/experiments/wil50_feasibility/fullrun/fullrun.par" -T 1 -i 1500 -u 48816 -v \
  -x "$PARENT" "$@" "$R/qaplib/wil50.dat" > "$U/$NAME.log" 2>&1
echo done > "$U/$NAME.done"
