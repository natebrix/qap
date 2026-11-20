#!/bin/bash
# arm17.sh NAME REBUILD PARENTFIX LINE: the depth-17 parent (its fix PARENTFIX, e.g. 8:49) from runs/rebuild/REBUILD/d16.out.state,
# forced to branch on LINE (QAP_FORCE_LINE=17:LINE; LINE "-": the params' own rule; PAR=file to use another params file), production's settings, its children bounded, depth 19 set aside (-d 18).
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-linetest
NAME=$1; RB=$2; PF=$3; LINE=$4; PAR=${PAR:-$B/infra/phase1/recipes/tai60b.par}; F=""; [ "$LINE" = "-" ] || F=17:$LINE; O=$D/runs/arm17/$NAME; mkdir -p $O
python3 -I $D/cutstate.py $D/runs/rebuild/$RB/d16.out.state $O/in.state 17 $PF
echo "./bnbdrv -p $PAR -r in.state -d 18 QAP_FORCE_LINE=$F $(date -u +%FT%TZ)" > $O/cmd.txt
cd $B/latest
set +e
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_FORCE_LINE=$F QAP_NODE_TRACE=$O/nodes.txt ./bnbdrv -p $PAR \
  --allow-settings-change -r $O/in.state -d 18 -o $O/out.state $B/qaplib/tai60b.dat > $O/log.txt 2> $O/err.txt
echo "exit $? $(date -u +%FT%TZ)" >> $O/cmd.txt
rm -f $O/in.state
