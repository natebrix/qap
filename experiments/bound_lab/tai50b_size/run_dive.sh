#!/bin/bash
# Knuth dives (bnbdrv -e) resumed from a state file of open nodes, one process on one core:
#   ./run_dive.sh SET SEED ARM IN.state CORE LEVEL SECONDS   -> runs/SET/sSEED/{est.out,dives.txt,dnn.txt,estimate.par}
# LEVEL = the depth of IN.state's nodes (estimate.par's start level: nothing above it is recomputed);
# stops at SECONDS of wall (a dive in flight is lost; every finished dive is in est.out and dives.txt).
H=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap; SET=$1; SEED=$2; ARM=$3; IN=$4; CORE=$5; LEVEL=$6; SECS=$7
D=$H/runs/$SET/s$SEED; mkdir -p $D $H/runs/store; cd $D
echo "$LEVEL 10000 1 1 1.5 1000" > estimate.par
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export QAP_EST_SEED=$SEED QAP_EST_HIST=$D/dives.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_WARM_DIR=$H/runs/store
CMD="$R/latest/bnbdrv -e -p $H/par/$ARM.par -u 458821518 -r $IN --allow-settings-change $R/qaplib/tai50b.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > cmd.txt
timeout -s TERM $SECS taskset -c $CORE $CMD > est.out 2> err.txt
echo "exit $?" >> cmd.txt; date -u +%FT%TZ >> cmd.txt
