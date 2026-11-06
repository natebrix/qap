#!/bin/bash
# Cold wil50 nodes of tuning/giveup/all32.txt, production flags with AMX (bundle
# D), one process at a time on core 1; speed/beta_sched/run_all32.sh with AMX and
# this directory's logs. NODES: line numbers of all32.txt (default: four
# fathoming nodes, 11 17 22 31). CONFIGS="name@env ..." as in run_warmset.sh.
# Logs: cold/<name>_<id>.log; summarize with speed/beta_sched/summarize.py's
# logic (summarize_cold.py here).
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CONFIGS=${CONFIGS:-"base@-"}
NODES=${NODES:-"11 17 22 31"}
CPU=${CPU:-1}
mkdir -p "$D/cold"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  for k in $NODES; do
    spec=$(sed -n "${k}p" "$D/../../tuning/giveup/all32.txt" | awk '{print $1}')
    f=$D/cold/${name}_$k.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=${AMX:-1} ${envs//+/ } \
      taskset -c "$CPU" "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -u 48816 -v ${EXTRA} \
      -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  done
done
