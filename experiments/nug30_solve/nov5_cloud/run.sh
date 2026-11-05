#!/bin/bash
# The Nov 4 record setting (laptop: 1,838 CPU s, 322 s wall, 3,844 nodes) on
# the cloud box: 4 workers, OpenBLAS, no AMX on this host. The checkpoint
# file ends with the final totals and per-level statistics (nug30_final.state).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
cd "$R"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_BETA=0.21 \
  QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) $(git rev-parse --short HEAD)" > "$D/host.txt"
python3 bench/parallel.py qaplib/nug30.dat --jobs 4 --ub 6124 --params PARAM/nug_deep_xbar_all.par \
  --extra "-x 1:7:1:200:1" --checkpoint "$D/nug30_final.state" --checkpoint-every 300 > "$D/run.log" 2>&1
