#!/bin/bash
# nug30 (m = 28) check of DNN_R32_NONCERT=1: eight depth-2 nodes of
# ../../../dnn_bound/predictor/nug30_depth2.csv (rows k = 3, 25, 75, 100, 125:
# fathoms; rows 0, 175: misses there) and the eigensolver study's 4:10,6:8,
# cold, with the nug30 record setting's DNN flags (PRSM 1.0/1.0, degree 3,
# -L 2, beta 0.21, refresh 120, one-pass CholQR), give-up after 200 as in the
# search, ub 6124. One process at a time on core $CORE, base then nc per node.
# PASS=2: the second pass, nodes in reverse order and nc first, logs
# nug30/base2_*.log and nug30/nc2_*.log.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CORE=${CORE:-0}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
NODES=${NODES:-"0:0,23:4 0:0,23:26 0:2,23:18 0:6,23:14 0:7,23:10 0:0,23:1 0:12,23:1 4:10,6:8"}
CONFS="base: nc:DNN_R32_NONCERT=1"; sfx=""
if [ "${PASS:-1}" = 2 ]; then
  NODES=$(echo $NODES | tr ' ' '\n' | tac | tr '\n' ' '); CONFS="nc:DNN_R32_NONCERT=1 base:"; sfx=2
fi
mkdir -p "$D/nug30"
for spec in $NODES; do
  tag=${spec//[:,]/_}
  for c in $CONFS; do
    name=${c%%:*}$sfx; envs=${c#*:}
    f=$D/nug30/${name}_$tag.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    env $envs taskset -c "$CORE" "$R/latest/dnntest" -T 1 -d 3 -L 2 -r 120 -B 0.21 -g 1.0 -G 1.0 \
      -H 200 -i 2000 -u 6124 -v -v -x "$spec" "$R/qaplib/nug30.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  done
done
