#!/bin/bash
# The root with the cut phase, bounded once and saved with its 50 children (state + parent store):
#   ./run_root.sh NAME CORE     -> runs/NAME/{log.txt,nodes.txt,dnn.txt,root.state,store/}
H=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap; NAME=$1; CORE=$2
D=$H/runs/$NAME; mkdir -p $D/store
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=$D/store
cd $R/latest
CMD="./bnbdrv -p $H/par/$NAME.par -u 458821518 -d 0 -o $D/root.state ../qaplib/tai50b.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/cmd.txt
taskset -c $CORE $CMD > $D/log.txt 2> $D/err.txt
echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
