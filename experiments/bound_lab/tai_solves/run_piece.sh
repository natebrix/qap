#!/bin/bash
# ./run_piece.sh NAME PAR INSTANCE UB CORE IN.state [EXTRA]: search every node of a state file to the end on one core
# -> runs/NAME/{log.txt,nodes.txt,dnn.txt,out.state (checkpoint every 600 s; no open nodes when done)}
H=$(cd "$(dirname "$0")" && pwd); NAME=$1; D=$H/runs/$NAME; mkdir -p $D $H/runs/store
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_WARM_DIR=$H/runs/store
cd /home/user/qap/latest
CMD="./bnbdrv -p $H/par/$2.par -u $4 -r $6 --allow-settings-change -o $D/out.state -k 600 ${7:-} ../qaplib/$3.dat"
{ echo "$CMD"; date -u +%FT%TZ; } >> $D/cmd.txt
taskset -c $5 $CMD >> $D/log.txt 2>> $D/err.txt; echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
