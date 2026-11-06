#!/bin/bash
# Default path unchanged: with none of the DNN_ADMM_* switches set, the binary
# built before the change (BASE) and after it (NEW) must print identical traces
# (all fields but the CPU columns) on the nug20 root (400 iterations, production
# flags, no cutoff) and on wil50 all32 node 19 (26:27,8:38,34:28,37:21, m = 46,
# fathomed at iteration 130). One process at a time on core CPU (1).
# Logs: ident/<name>_{base,new}.log; prints OK or DIFF per case.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BASE=${BASE:?path of the pre-change dnntest}
NEW=${NEW:-$R/latest/dnntest}
CPU=${CPU:-1}
mkdir -p "$D/ident"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=${AMX:-1}
unset DNN_ADMM_CROSS DNN_ADMM_WARM_BETA DNN_ADMM_WARM_DAMP
strip() { sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+//; s/ ms_per_iter=[0-9.]+//; /^profile|^eig split|^step split|^form R|^refresh reasons/d' "$1"; }
run() {  # name flags...
  local name=$1; shift
  for v in base new; do
    local bin=$BASE; [ $v = new ] && bin=$NEW
    taskset -c "$CPU" "$bin" "$@" > "$D/ident/${name}_$v.log" 2>&1
  done
  if cmp -s <(strip "$D/ident/${name}_base.log") <(strip "$D/ident/${name}_new.log"); then echo "$name: OK ($(grep -c '^it ' "$D/ident/${name}_new.log") trace lines)"; else echo "$name: DIFF"; fi
}
run nug20_root -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -i 400 -v "$R/qaplib/nug20.dat"
run wil50_n19 -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -u 48816 -v -x 26:27,8:38,34:28,37:21 "$R/qaplib/wil50.dat"
# a warm child with the production beta setting (the warm path through BetaSched)
run wil50_warm -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 1500 -u 48816 -v -x 20:12,40:28,17:25 -K 43:18 -Q b "$R/qaplib/wil50.dat"
