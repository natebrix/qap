#!/bin/bash
# ./run.sh NAME BIN ARGS...: one bnbdrv process (single-threaded) in runs/NAME with node and DNN traces
# BIN: orig (cognac 0e6c2ecf as built) or hook (with the QAP_MIRROR hook, mirror.patch). Env QAP_MIRROR passes through.
H=$(cd "$(dirname "$0")" && pwd); NAME=$1; BIN=$2; shift 2; D=$H/runs/$NAME; mkdir -p $D
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=${STORE:-$H/runs/store}
mkdir -p $QAP_DNN_WARM_DIR
cd /home/user/wt-mirror/latest
B=./bnbdrv.$BIN
{ echo "QAP_MIRROR=${QAP_MIRROR:-} $B $*"; date -u +%FT%TZ; } > $D/cmd.txt
$B "$@" > $D/log.txt 2> $D/err.txt; echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
