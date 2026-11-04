#!/bin/bash
# E6: end-to-end A/B in the search on the sub-subtree of [9:45] (node 10 +
# 9:45, depth 4): its cold first solve, its cut phase (a miss, 63.8 short in
# pool_search), then its 46 children on facility 19's row, which inherit the
# pool and the phase's state (the pooled class, about a quarter of the
# node-10 subtree's DNN CPU). Exactly run_pool_search.sh's configuration,
# serial (one bnbdrv, no checkpoints, so every warm start is kept).
# Arm plain: as production. Arm polish: plus DNN_CERTX=32 DNN_CERTX_USE=1
# (the reduce-pass skew polish, 12 passes + 1 exact, used as the certificate
# at checks within W of the cutoff; lean: the standard lambda bound).
# Usage: run_e6.sh <core> plain|polish [W]
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; arm=$2; W=${3:-10}
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
name=e6_n9_45_$arm; [ "$arm" = polish ] && name=${name}_w$W
export QAP_DNN_TRACE=$D/logs/$name.trace
if [ "$arm" = polish ]; then
  export DNN_CERTX=32 DNN_CERTX_USE=1 DNN_CERTX_QUIET=1 DNN_CERTX_EVERY=1000000 DNN_CERTX_WINDOW=$W
fi
f=$D/logs/$name.log
grep -q "Total Time" "$f" 2>/dev/null && exit 0
rm -f "$QAP_DNN_TRACE"
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) core $core $(cd "$R" && git rev-parse --short HEAD)" > "$f.tmp"
(cd "$R/latest" && taskset -c "$core" ./bnbdrv -p PARAM/nug_deep_xbar9_all.par -F 39:27,30:48,25:49,9:45 \
   -x 0:4:1:300:1 -u 48816 "$R/qaplib/wil50.dat") >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
