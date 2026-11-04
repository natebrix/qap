#!/bin/bash
# Pool inheritance in the search, measured (README.md, predeclared). The
# stacked arm of stronger_bounds/search/run_ab6.sh (node-10 subtree, rule 9 at
# every node through depth 6, cut phase zone 0.003, quiet 1,500, stall 3, EXT
# 20, bundle D, QAP_DNN_WARM_Y=2; 322 nodes, 32,740 CPU-s,
# wil50n10ext_cut_rule9_stacked.log) with one change: QAP_DNN_POOL=1 (a node
# whose cut phase beat its uncut bound hands its pool and that solve's state
# to its children; latest/dnncut.h). QAP_DNN_POOL_CHECK=20 guards the
# assembled terms; QAP_DNN_TRACE logs a line per DNN node (trace.txt).
# 4 workers, checkpoint every 30 min (the pools survive a checkpoint, the DNN
# states do not, as in the stacked arm). The host is logged at every
# (re)start (hosts.txt). Resumable (Total Time line = done; else the .ckpt).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2
export QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20 QAP_DNN_TRACE=$D/trace.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
name=wil50n10ext_cut_rule9_stacked_pool; f=$D/$name.log
grep -q "Total Time" "$f" 2>/dev/null && exit 0
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') amx_cpus $(grep -c amx_bf16 /proc/cpuinfo) $(cd $R && git rev-parse --short HEAD)" >> $D/hosts.txt
resume=(); [ -s "$D/$name.ckpt" ] && resume=(--resume "$D/$name.ckpt")
(cd "$R" && python3 bench/parallel.py qaplib/wil50.dat --jobs 4 --ub 48816 \
   --params PARAM/nug_deep_xbar9_all.par --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" \
   --checkpoint "$D/$name.ckpt" --checkpoint-every 1800 "${resume[@]}") >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
