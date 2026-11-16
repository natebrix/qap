#!/bin/bash
# ./run_dive.sh SET SEED PAR INSTANCE UB CORE IN.state LEVEL TRIALS SECONDS: Knuth dives (bnbdrv -e -r, est_resume.patch)
# from a state file's nodes (LEVEL = their depth), TRIALS dives or SECONDS of wall -> runs/SET/sSEED/
H=$(cd "$(dirname "$0")" && pwd); D=$H/runs/$1/s$2; mkdir -p $D $H/runs/store; cd $D
echo "$8 $9 1 1 1.5 1000" > estimate.par
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_EST_SEED=$2 QAP_EST_HIST=$D/dives.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_WARM_DIR=$H/runs/store
CMD="/home/user/qap/latest/bnbdrv -e -p $H/par/$3.par -u $5 -r $7 --allow-settings-change /home/user/qap/qaplib/$4.dat"
{ echo "$CMD"; date -u +%FT%TZ; } > cmd.txt
timeout -s TERM ${10} taskset -c $6 $CMD > est.out 2> err.txt; echo "exit $?" >> cmd.txt; date -u +%FT%TZ >> cmd.txt
