#!/bin/bash
# dumpnode.sh NAME SRC.state DEPTH FIX: what rule 9 saw at a node of the c47 descent (README, "Rule 9's failures"):
# the node, cut from SRC (its parent's state, children set aside), bounded and branched under tai60b.par (rule 9 below
# depth 6) with its children set aside (-d DEPTH), dumping its QPB U, bound, cutoff and DNN marginals x-bar
# (QAP_LINE_SCORES=DEPTH), every line's lookahead score (QAP_LINE_SCORES_LA=1), every line's marginal score ranked
# (QAP_XBAR_DEBUG=2) and rule 9's arbitration (xbar9). Output in runs/dump/NAME/err.txt.
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-linetest
NAME=$1; SRC=$2; DEP=$3; FIX=$4; O=$D/runs/dump/$NAME; mkdir -p $O
python3 -I $D/cutstate.py $SRC $O/in.state $DEP $FIX
cd $B/latest; set +e
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_LINE_SCORES=$DEP QAP_LINE_SCORES_LA=1 QAP_XBAR_DEBUG=2 QAP_NODE_TRACE=$O/nodes.txt \
  ./bnbdrv -p $B/infra/phase1/recipes/tai60b.par --allow-settings-change -r $O/in.state -d $DEP $B/qaplib/tai60b.dat > $O/log.txt 2> $O/err.txt
echo "exit $?" > $O/done.txt; rm -f $O/in.state
