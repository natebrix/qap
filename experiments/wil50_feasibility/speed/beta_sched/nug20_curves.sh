#!/bin/bash
# nug20 root, cold, production flags (degree 3, -L 2, PRSM, -r 40, one CholQR
# pass), 400 iterations, no cutoff and no give-up, one per configuration on
# core 1. CONFIGS="name@env ..." as in run_all32.sh. Logs: nug20/<name>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CONFIGS=${CONFIGS:-"base@-"}
INST=${INST:-nug20}
ITERS=${ITERS:-400}
FLAGS=${FLAGS:-"-d 3 -L 2 -g 1.0 -G 1.0 -r 40"}
mkdir -p "$D/$INST"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 ${envs//+/ } \
    taskset -c 1 "$BIN" -T 1 $FLAGS -H 0 -i "$ITERS" -v $EXTRA "$R/qaplib/$INST.dat" > "$D/$INST/$name.log" 2>&1
done
