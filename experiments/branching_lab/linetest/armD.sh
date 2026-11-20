#!/bin/bash
# armD.sh NAME SRC.state DEPTH FIX PAR: the depth-DEPTH node with fix FIX from SRC (a state holding one parent's
# children set aside), bounded under PAR with the params' own line rule, its children bounded, the next level set aside
# (-d DEPTH+1). Runs the latest_dbg build (latest/ + maxf_debug.patch: QAP_MAXF_DEBUG=1 prints rule 2's counts; the
# search is unchanged). For the rule 2 / rule 9 scan (README).
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-linetest
NAME=$1; SRC=$2; DEP=$3; FIX=$4; PAR=$5; O=$D/runs/scan/$NAME; mkdir -p $O
python3 -I $D/cutstate.py $SRC $O/in.state $DEP $FIX
echo "./bnbdrv -p $PAR -r in.state -d $((DEP+1)) (latest_dbg, QAP_MAXF_DEBUG=1) $(date -u +%FT%TZ)" > $O/cmd.txt
cd $B/latest_dbg; set +e
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_MAXF_DEBUG=1 QAP_NODE_TRACE=$O/nodes.txt ./bnbdrv -p $PAR --allow-settings-change -r $O/in.state \
  -d $((DEP+1)) -o $O/out.state $B/qaplib/tai60b.dat > $O/log.txt 2> $O/err.txt
echo "exit $? $(date -u +%FT%TZ)" >> $O/cmd.txt; rm -f $O/in.state
