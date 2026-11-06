#!/bin/bash
# Scan the AMX kernel's runtime knobs (DNN_CORE_AMX_IB/KB/NT/ORDER, latest/amx_gemm.c)
# with ../../amx/test_amx_gemm.c at K = 2117, kb = 610 (the profile's mean block), on
# CPU 3, one process at a time. Every variant is interleaved with the default (base),
# so a drift of the host's load hits both; 3 runs each, the kernel's own timer
# (ms per product: Y split + tile loop) and the OpenBLAS sgemm reference.
D=$(cd "$(dirname "$0")" && pwd)
CPU=${CPU:-3}
K=${K:-2117}; KB=${KBLOCK:-610}
OUT=${OUT:-$D/logs/scan_amx.txt}
echo "# $(date -u +%FT%TZ) K $K kb $KB load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$OUT"
run() { # name env...
  local name=$1; shift
  for r in 1 2 3; do
    local line
    line=$(env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$@" taskset -c $CPU "$D/test_amx" $K $KB | tr '\n' ' ')
    echo "$name run$r: $line" >> "$OUT"
  done
}
VARIANTS=${VARIANTS:-"nt0:DNN_CORE_AMX_NT=0 ord1:DNN_CORE_AMX_ORDER=1 nt0ord1:DNN_CORE_AMX_NT=0+DNN_CORE_AMX_ORDER=1 ib64:DNN_CORE_AMX_IB=64 ib256:DNN_CORE_AMX_IB=256 kb4:DNN_CORE_AMX_KB=4 kb16:DNN_CORE_AMX_KB=16 kb67:DNN_CORE_AMX_KB=67 ib64kb16:DNN_CORE_AMX_IB=64+DNN_CORE_AMX_KB=16 ib256kb4:DNN_CORE_AMX_IB=256+DNN_CORE_AMX_KB=4"}
run base
for v in $VARIANTS; do
  name=${v%%:*}; envs=${v#*:}
  run "$name" ${envs//+/ }
  run base
done
echo "# done $(date -u +%T)" >> "$OUT"
