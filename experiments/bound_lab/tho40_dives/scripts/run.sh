#!/bin/bash
# One estimator process: scripts/run.sh NAME PAR STRATUM NDIVES SEED UB START RULE [FULL_DEPTH] [PATH_PERM]
#   STRATUM  the depth-1 node 26:c the process starts from (-F), bounded once, exactly
#   START    estimate.par's start level: the depth at which dives start (2: from the depth-2 children)
#   RULE     qpb: children drawn by their own QPB bound (QAP_EST_CHILD_QPB=1); pred: by the predicted bound
#   FULL_DEPTH  (patched binary) from this depth down a dive takes every child (the subtree is enumerated); - for none
#   PATH_PERM   (patched binary) a dive follows this permutation instead of drawing
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
BIN=${DIVE_BIN:-/home/user/dvbuild/bnbdrv}; [ "$2" = hand ] && BIN=/home/user/dvbuild3/bnbdrv
NAME=$1 PAR=$2 STRATUM=$3 ND=$4 SEED=$5 UB=$6 START=$7 RULE=$8 FULL=${9:--} PPERM=${10:-}
D=$HERE/runs/$NAME; mkdir -p "$D"; cd "$D" || exit 1
echo "$START $ND 1 100000000 10 1000" > estimate.par
export OPENBLAS_NUM_THREADS=1 QAP_EST_SEED=$SEED QAP_EST_HIST=$D/dives.txt QAP_DNN_TRACE=$D/trace.txt DIVE_TRACE=$D/trace.txt
[ "$RULE" = qpb ] && export QAP_EST_CHILD_QPB=1
[ "$FULL" != - ] && export DIVE_FULL_DEPTH=$FULL
[ -n "$PPERM" ] && export DIVE_PATH="$PPERM"
CMD="$BIN -e -p $HERE/par/$PAR.par -s $ROOT/latest/SYM/tho40t.sym -u $UB -F $STRATUM $ROOT/qaplib/tho40t.dat"
{ echo "start $(date -u +%FT%TZ) load $(cut -d' ' -f1-3 /proc/loadavg)"; echo "env QAP_EST_SEED=$SEED QAP_EST_CHILD_QPB=${QAP_EST_CHILD_QPB:-0} DIVE_FULL_DEPTH=${DIVE_FULL_DEPTH:-} DIVE_PATH=${DIVE_PATH:-}"; echo "estimate.par: $(cat estimate.par)"; echo "$CMD"; } > cmd.txt
${RUN_PREFIX:-} $CMD > est.out 2> err.txt
echo "exit $? $(date -u +%FT%TZ) load $(cut -d' ' -f1-3 /proc/loadavg)" >> cmd.txt
