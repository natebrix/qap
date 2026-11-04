#!/bin/bash
# Smoke test of DNN_CERTX_USE=1 (the polished value used for fathoming): the
# same node with mode 0 (plain), 2 (exact-pass polish), 32 (reduce-pass polish)
# and 34, evaluated within <window> of the cutoff only.
# Usage: use_smoke.sh <dnntest> <instance> <ub> <window> [fixes]
exe=$1; inst=$2; ub=$3; win=$4; fx=${5:--}
for mode in 0 2 32 34; do
  DNN_CERTX=$mode DNN_CERTX_USE=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=$win OMP_NUM_THREADS=1 \
    "$exe" -T 1 -g 1.0 -G 1.0 -H 0 -i 1500 -u "$ub" -x "$fx" "$inst" | grep -E "^RESULT" | sed "s/^/mode $mode: /" | cut -c1-170
done
