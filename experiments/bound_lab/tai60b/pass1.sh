#!/bin/bash
# pass1.sh TAG INSTATE DEPTH K SEED [F:L ...]: the first (pooled, warm) solves of K random open nodes of INSTATE
# (plus the forced ones), no phase (par/kid0.par: cut.zone = 0): which fathom at once, how short are the rest.
# (after regions/size50/pass1.sh). SKIP/OF in the environment: take every OF-th sampled node starting at SKIP,
# so that one sample can be split over processes (TAG then gets the suffix _SKIP).
D=$(cd $(dirname $0); pwd); I=${INST:-tai60b}
TAG=$1; IN=$2; DEPTH=$3; K=$4; SEED=$5; shift 5
mkdir -p $D/runs/$TAG; rm -f $D/runs/$TAG/nodes.txt $D/runs/$TAG/dnn.txt
python3 -I $D/cutset.py $IN $D/runs/$TAG/in.state $K $SEED "$@" > $D/runs/$TAG/sample.txt || exit 1
cd /home/user/wt-regions/latest
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_NODE_TRACE=$D/runs/$TAG/nodes.txt QAP_DNN_TRACE=$D/runs/$TAG/dnn.txt
CMD="./bnbdrv -p $D/par/${PAR0:-kid0}.par --allow-settings-change -r $D/runs/$TAG/in.state -d $DEPTH -o $D/runs/$TAG/out.state ../qaplib/$I.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/runs/$TAG/cmd.txt
$CMD > $D/runs/$TAG/log.txt 2> $D/runs/$TAG/err.txt
{ echo "exit $?"; date -u +%FT%TZ; uptime; } >> $D/runs/$TAG/cmd.txt
rm -f $D/runs/$TAG/in.state $D/runs/$TAG/out.state
