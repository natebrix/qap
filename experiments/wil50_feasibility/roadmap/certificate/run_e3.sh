#!/bin/bash
# E3: end-to-end A/B of the polished certificate used for fathoming
# (DNN_CERTX_USE=1, lean path: the standard certificate's lambda, only the LAP
# passes added), evaluated only once the running best is within W of the
# cutoff. Arm A: plain dnntest; arm B: DNN_CERTX=<mode> USE=1. Same flags as
# E1 (cold depth-4 children of node 10, production solver settings).
# Usage: run_e3.sh <core> <mode> <window> <child> [<child> ...]
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; mode=$2; win=$3; shift 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
run() {   # name, env...
  local f=$D/logs/$1.log; shift
  grep -q '^RESULT' "$f" 2>/dev/null && return
  echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core fixes $fx" > "$f.tmp"
  env "$@" taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 \
    -x "$fx" "$R/qaplib/wil50.dat" >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
# a child given with ':' ("49:45,0:12") is a list of extra fixes under node 10
# (the Nov 4 run passed it through the X:45 template; dnntest read "0:12:45"
# as 0:12, so those logs are for the intended node)
for x in "$@"; do
  case $x in *:*) fx="39:27,30:48,25:49,$x";; *) fx="39:27,30:48,25:49,$x:45";; esac
  n=$(echo "$x" | tr ':,' '-_')
  run e3_c${n}_plain DNN_CERTX=0
  run e3_c${n}_m${mode}_w${win} DNN_CERTX=$mode DNN_CERTX_USE=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=$win
done
