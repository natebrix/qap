#!/bin/bash
# xdump.sh NODE: re-bound one node whose lines were bounded here, under the same settings, with its children set
# aside, and dump its DNN marginals x-bar (QAP_XBAR_DUMP, in the copy /home/user/wt-cognac/latest_x: lsla.patch and xbar_dump.patch;
# plus the dump). NODE: d7 (10:39, as run.sh strong), d8 (14:54, as at8), d9 (24:2, as at9), r50 (tai50b's root,
# as root50), r35 (tai35b's root, as t2). Output runs/xdump/NODE/{xbar,nodes,log,err}.txt.
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd); B=/home/user/wt-cognac; X=$B/latest_x
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
N=$1; O=$D/runs/xdump/$N; mkdir -p $O; rm -f $O/{xbar,nodes,log,err,done}.txt $O/in.state
T60=$B/qaplib/tai60b.dat
case $N in
  d7) python3 -I $D/cutstate.py $D/runs/parent/out.state $O/in.state 7 10:39; DEP=7
      A=(-p $D/strong.par --allow-settings-change -r $O/in.state -d 7 $T60) ;;
  d8) python3 -I $D/cutstate.py $D/runs/set8_c54/out.state $O/in.state 8 14:54; DEP=8
      A=(-p $D/strong.par --allow-settings-change -r $O/in.state -d 8 $T60) ;;
  d9) python3 -I $D/cutstate.py $D/runs/set9_c2/out.state $O/in.state 9 24:2; DEP=9
      A=(-p $D/strong.par --allow-settings-change -r $O/in.state -d 9 $T60) ;;
  r50) DEP=0; A=(-p $B/infra/phase1/recipes/tai50b.par -u 458821518 -d 0 $B/qaplib/tai50b.dat) ;;
  r35) DEP=0; A=(-p $D/tai35b.par -u 283315446 -d 0 $B/qaplib/tai35b.dat) ;;
  *) echo "NODE: d7 d8 d9 r50 r35"; exit 2 ;;
esac
cd $X; set +e
QAP_ALLOW_UNKNOWN=1 QAP_XBAR_DUMP=$DEP:$O/xbar.txt QAP_NODE_TRACE=$O/nodes.txt ./bnbdrv "${A[@]}" > $O/log.txt 2> $O/err.txt
echo "exit $? $(date -u +%FT%TZ)" > $O/done.txt; rm -f $O/in.state
