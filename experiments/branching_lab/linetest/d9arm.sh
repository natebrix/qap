#!/bin/bash
# d9arm.sh NAME SRC.state FIX PAR: the depth-8 node with fix FIX from SRC (a strong7 state: 10:39 branched, its depth-8
# children set aside), bounded under PAR with the params' own line rule, its depth-9 children bounded, depth 10 set
# aside (-d 9). For the 2x2 at depth 9 (README).
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-linetest
NAME=$1; SRC=$2; FIX=$3; PAR=$4; O=$D/runs/d9x/$NAME; mkdir -p $O
python3 -I $D/splitstate.py $SRC $O/in.state 8 $FIX 2> $O/split.txt
echo "./bnbdrv -p $PAR -r in.state -d 9 $(date -u +%FT%TZ)" > $O/cmd.txt
cd $B/latest; set +e
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_NODE_TRACE=$O/nodes.txt ./bnbdrv -p $PAR --allow-settings-change -r $O/in.state -d 9 \
  -o $O/out.state $B/qaplib/tai60b.dat > $O/log.txt 2> $O/err.txt
echo "exit $? $(date -u +%FT%TZ)" >> $O/cmd.txt; rm -f $O/in.state
