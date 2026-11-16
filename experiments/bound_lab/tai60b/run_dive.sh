#!/bin/bash
# run_dive.sh SET SEED IN.state LEVEL TRIALS [SECONDS]: Knuth dives (bnbdrv -e -r, after ../tai50b_size/run_dive.sh)
# from the open nodes of IN.state (all at depth LEVEL), one process on one core, under par/kidW.par with the
# incumbent at best known + 1. Stops after TRIALS finished dives (or SECONDS of wall, losing the dive in flight).
# -> runs/SET/sSEED/{est.out,dives.txt,dnn.txt,estimate.par,cmd.txt}
H=$(cd "$(dirname "$0")" && pwd); SET=$1; SEED=$2; IN=$3; LEVEL=$4; TRIALS=$5; SECS=${6:-86400}
D=$H/runs/$SET/s$SEED; mkdir -p $D $H/runs/store; cd $D; rm -f dives.txt dnn.txt
echo "$LEVEL $TRIALS 1 1 1.5 1000" > estimate.par
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export QAP_EST_SEED=$SEED QAP_EST_HIST=$D/dives.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=$H/runs/store
CMD="/home/user/wt-regions/latest/bnbdrv -e -p $H/par/${PAR:-kidW}.par -u 608215055 -r $IN --allow-settings-change /home/user/qap/qaplib/tai60b.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > cmd.txt
timeout -s TERM $SECS $CMD > est.out 2> err.txt
echo "exit $?" >> cmd.txt; date -u +%FT%TZ >> cmd.txt
