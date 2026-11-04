#!/bin/bash
# The stacked arm WITHOUT pool inheritance through the coordinator, on the
# same host as run_coord_search.sh (the pool arm), to measure inheritance
# like for like: the parallel.py comparison (pool_search, 1.37x) used a
# stacked arm run on Nov 3-4 hosts of unknown speed, and the host after the
# Nov 4 17:49 restart is 1.23x faster per iteration than the one before it.
# Same settings as run_coord_search.sh except QAP_DNN_POOL (off).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
B=local:$D/backend_stacked; RUN=n10stacked
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$R"
grep -q "coordinator: 4 slots" "$D/coordinator_stacked.log" 2>/dev/null && exit 0
echo "$(date -u +%FT%TZ) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') amx_cpus $(grep -c amx_bf16 /proc/cpuinfo) $(git rev-parse --short HEAD)" >> "$D/hosts_stacked.txt"
if [ ! -f "$D/backend_stacked/ledger.db" ]; then
  python3 -m coord.qapctl init $RUN --backend $B --instance qaplib/wil50.dat --params PARAM/nug_deep_xbar9_all.par \
    --ub 48816 --slots 4 --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --trace \
    --budget 3600 --target-unit-cpu 1800 --steal-after 10 --lease 600 --heartbeat 60 --ckpt-upload 300 \
    --checkpoint-every 60 --status-every 60 >> "$D/coordinator_stacked.log" 2>&1 || exit 1
else
  python3 -m coord.qapctl requeue-leased $RUN --backend $B >> "$D/coordinator_stacked.log" 2>&1
fi
python3 -m coord.qapctl worker $RUN --backend $B --slots 4 --workdir "$D/work_stacked" >> "$D/worker_stacked.log" 2>&1 &
python3 -m coord.qapctl coordinator $RUN --backend $B >> "$D/coordinator_stacked.log" 2>&1
wait
