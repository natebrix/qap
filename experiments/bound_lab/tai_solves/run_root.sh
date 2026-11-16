#!/bin/bash
# ./run_root.sh NAME PAR INSTANCE UB CORE: the root bounded with its phase and saved with its children -> runs/NAME/
H=$(cd "$(dirname "$0")" && pwd); NAME=$1; D=$H/runs/$NAME; mkdir -p $D/store
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=$D/store
cd /home/user/qap/latest
CMD="./bnbdrv -p $H/par/$2.par -u $4 -d 0 -o $D/root.state ../qaplib/$3.dat"
{ echo "$CMD"; date -u +%FT%TZ; } > $D/cmd.txt
taskset -c $5 $CMD > $D/log.txt 2> $D/err.txt; echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
