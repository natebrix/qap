#!/bin/bash
# A/B of the ADMM penalty schedules (DNN_BETA_SCHED, latest/dnn.C) on the 32
# wil50 nodes of tuning/giveup/all32.txt, cold, in the search's production
# mode: degree 3, tight lambda (-L 2), PRSM (-g 1.0 -G 1.0), a full
# eigendecomposition every 40 iterations (-r 40), give-up after 300 (-H 300),
# -i 2500 -u 48816, one CholQR pass (DNN_CHOLQR_PASSES=1), one BLAS thread,
# ONE run at a time on core 1 (taskset -c 1; the other cores are in use).
#
# Configurations: CONFIGS="name@env ..." (env: space-free assignments joined
# by '+', or '-' for none). NODES: the line numbers of all32.txt to run, in
# order (default 1..32). BIN: the dnntest binary (default latest/dnntest).
# Resumable: a node whose log has a RESULT line is skipped.
# Logs: runs/<name>_<id>.log, id = line number in all32.txt
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CONFIGS=${CONFIGS:-"base@-"}
NODES=${NODES:-$(seq 1 32)}
mkdir -p "$D/runs"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  for k in $NODES; do
    spec=$(sed -n "${k}p" "$D/../../tuning/giveup/all32.txt" | awk '{print $1}')
    f=$D/runs/${name}_$k.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 ${envs//+/ } \
      taskset -c 1 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -u 48816 -v \
      -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  done
done
