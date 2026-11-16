#!/bin/bash
# Resume a state file under par/ARM.par on one core, bounding its nodes down to depth D and saving the open nodes below:
#   ./run_state.sh NAME ARM IN.state CORE D [EXTRA]   -> runs/NAME/{log.txt,nodes.txt,dnn.txt,out.state}; parent states in runs/store
H=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap; NAME=$1; ARM=$2; IN=$3; CORE=$4; DEP=$5; EXTRA=${6:-}
D=$H/runs/$NAME; mkdir -p $D $H/runs/store
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=$H/runs/store
cd $R/latest
CMD="./bnbdrv -p $H/par/$ARM.par -u 458821518 -r $IN --allow-settings-change -d $DEP -o $D/out.state $EXTRA ../qaplib/tai50b.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/cmd.txt
taskset -c $CORE $CMD > $D/log.txt 2> $D/err.txt
echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
