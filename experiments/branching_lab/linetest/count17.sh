#!/bin/bash
# count17.sh: rule 2's count (QAP_MAXF_DEBUG) at the three depth-17 parents of the line test: each parent bounded and
# branched under rule2.par with the latest_dbg build, its children set aside (-d 17), one at a time.
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-linetest; RB=$D/runs/rebuild
for x in "p849 a849 8:49" "p857 a849 8:57" "p2231 path 22:31"; do set -- $x; O=$D/runs/count17/$1; mkdir -p $O
  python3 -I $D/cutstate.py $RB/$2/d16.out.state $O/in.state 17 $3
  (cd $B/latest_dbg && OPENBLAS_NUM_THREADS=1 QAP_MAXF_DEBUG=1 QAP_NODE_TRACE=$O/nodes.txt ./bnbdrv -p $D/rule2.par --allow-settings-change \
     -r $O/in.state -d 17 $B/qaplib/tai60b.dat > $O/log.txt 2> $O/err.txt); echo "exit $?" > $O/done.txt; rm -f $O/in.state
done
