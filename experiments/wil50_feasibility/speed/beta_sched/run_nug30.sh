#!/bin/bash
# nug30 depth-2 nodes (nug30_nodes.txt: seven from dnn_bound/predictor/
# nug30_depth2.csv, two misses and five fathoms of 100-320 iterations there,
# and 4:10,6:8 of experiments/eigensolver/README.md), cold, the nug30
# record's DNN setting: degree 3, -L 2, PRSM, QAP_DNN_BETA=0.21 (-B 0.21),
# QAP_DNN_REFRESH=120 (-r 120), one CholQR pass; -H 300 -i 1000, incumbent
# 6124. One process on core 1 per configuration. CONFIGS="name@env ..." as
# in run_all32.sh. NODEFILE (nug30_nodes.txt) and SUFFIX select another node
# list: nug30_nodes_big.txt holds the 28 depth-2 nodes k = 4 mod 8 of that csv.
# Logs: nug30/<name><SUFFIX>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CONFIGS=${CONFIGS:-"base@-"}
mkdir -p "$D/nug30"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 ${envs//+/ } \
    taskset -c 1 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -B 0.21 -r 120 -H 300 -i 1000 -u 6124 -v \
    -N "$D/${NODEFILE:-nug30_nodes.txt}" "$R/qaplib/nug30.dat" > "$D/nug30/$name$SUFFIX.log" 2>&1
done
