#!/bin/bash
# The whole search, level by level, on 4 cores:  ./bfs.sh TAG PAR INSTANCE UB START.state FROM UNTIL [CHUNKS]
# At depth D = FROM..UNTIL-1 the open nodes (all at depth D) are dealt into CHUNKS state files and each is
# bounded with -d D by one bnbdrv process (4 at a time), whose -o file holds the children that stay open; the
# files are merged (splitstate.py -merge) for depth D+1. At depth UNTIL the chunks are searched to the end.
# Every node starts warm from its parent with its pool, through the state files and the parent store.
# runs/TAG/dD_cK/{log.txt,nodes.txt,dnn.txt,cmd.txt}; runs/TAG/open_D.state; runs/TAG/progress.txt
H=$(cd "$(dirname "$0")" && pwd); cd $H; TAG=$1; PAR=$2; INST=$3; UB=$4; CUR=$5; FROM=$6; UNTIL=$7; K=${8:-8}
R=runs/$TAG; mkdir -p $R $R/store
# the start file names its parents' DNN states relative to its own directory when its store is below it: link them here
[ -d "$(dirname $CUR)/store" ] && ln -sf $H/$(dirname $CUR)/store/*.dnn $R/store/ 2>/dev/null
one() {  # NAME IN EXTRA
  D=$H/$R/$1; mkdir -p $D
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_WARM_DIR=$H/$R/store
  cd /home/user/qap/latest
  CMD="./bnbdrv -p $H/par/$PAR.par -u $UB -r $2 --allow-settings-change -o $D/out.state $3 ../qaplib/$INST.dat"
  { echo "$CMD"; date -u +%FT%TZ; } > $D/cmd.txt
  $CMD > $D/log.txt 2> $D/err.txt; echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
}
export -f one; export H R PAR INST UB
for d in $(seq $FROM $UNTIL); do
  n=$(grep -c "^node " $CUR); echo "$(date -u +%T) depth $d: $n open nodes" >> $R/progress.txt
  [ "$n" -eq 0 ] && break
  k=$K; [ "$n" -lt "$k" ] && k=$n
  python3 splitstate.py $CUR $R/in_${d}_ $k > /dev/null
  extra="-d $d"; [ "$d" -eq "$UNTIL" ] && extra=""
  seq 0 $((k-1)) | xargs -P ${NPROC:-4} -I{} bash -c "one d${d}_c{} $H/$R/in_${d}_{}.state '$extra'"
  [ "$d" -eq "$UNTIL" ] && break
  python3 splitstate.py -merge $R/open_$((d+1)).state $(for i in $(seq 0 $((k-1))); do echo $R/d${d}_c$i/out.state; done) > /dev/null
  CUR=$R/open_$((d+1)).state
done
echo "$(date -u +%T) done" >> $R/progress.txt
