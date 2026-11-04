#!/bin/bash
# The node-10 subtree through the coordinator (coord/, plan.md phase 1):
# exactly pool_search/run_pool_search.sh's configuration (stacked arm + pool
# inheritance; 322 nodes, 23,905 CPU-s, wall 7,342 s with bench/parallel.py,
# 4 workers, checkpoint barrier every 30 minutes), with the coordinator and
# one 4-slot worker instead of parallel.py. Units: CPU budget 3,600 s before a
# unit returns its open nodes, new units of about 1,800 CPU-s, 2 x slots kept
# queued, an idle slot may ask a unit that has run 10 s to return. Lease 600 s,
# heartbeat 60 s, the solver's checkpoint (-k 60) uploaded every 300 s.
# Resumable: after a container restart run it again (the coordinator replays
# from the ledger, the dead worker's units are re-queued and resume from their
# checkpoints). The host is logged at every (re)start (hosts.txt).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
B=local:$D/backend; RUN=n10pool
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$R"
grep -q "coordinator: 4 slots" "$D/coordinator.log" 2>/dev/null && exit 0
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') amx_cpus $(grep -c amx_bf16 /proc/cpuinfo) $(git rev-parse --short HEAD)" >> "$D/hosts.txt"
if [ ! -f "$D/backend/ledger.db" ]; then
  python3 -m coord.qapctl init $RUN --backend $B --instance qaplib/wil50.dat --params PARAM/nug_deep_xbar9_all.par \
    --ub 48816 --slots 4 --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --trace \
    --budget 3600 --target-unit-cpu 1800 --steal-after 10 --lease 600 --heartbeat 60 --ckpt-upload 300 \
    --checkpoint-every 60 --status-every 60 >> "$D/coordinator.log" 2>&1 || exit 1
else
  python3 -m coord.qapctl requeue-leased $RUN --backend $B >> "$D/coordinator.log" 2>&1
fi
python3 -m coord.qapctl worker $RUN --backend $B --slots 4 --workdir "$D/work" >> "$D/worker.log" 2>&1 &
python3 -m coord.qapctl coordinator $RUN --backend $B >> "$D/coordinator.log" 2>&1
wait
