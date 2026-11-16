#!/bin/bash
# tai50b level by level on 4 cores:  ./bfs.sh START.state FROM UNTIL
# At depth D = FROM..UNTIL the open nodes (all at depth D) are cut into up to 48 contiguous state files (one
# parent's children together), each bounded with -d D by one bnbdrv process (4 at a time); the children left open
# are merged into runs/open_{D+1}.state. After each level: open/open_{D+1}.txt (the open nodes' fixes, committed;
# the state files are not), the parents' DNN states no longer referenced are deleted, levels.py rewrites the README
# table and checks the stop conditions (a non-zero exit stops the run), and the level is committed and pushed.
# Rerunning an interrupted level (./bfs.sh runs/open_D.state D 48) redoes only its unfinished chunks.
# A file named "stop" in this directory stops the run before the next level.
# Binary: $BNBDIR/bnbdrv = latest + regions.patch + est_resume.patch + pool_retry.patch (both steps).
H=$(cd "$(dirname "$0")" && pwd); cd $H; CUR=$1; FROM=$2; UNTIL=$3
export BNBDIR=${BNBDIR:-/home/user/qap_retry2/latest} UB=458821518 PAR=A_wc INST=tai50b H
mkdir -p runs/store
[ -d "$(dirname $CUR)/store" ] && [ "$(dirname $CUR)" != runs ] && ln -sf $(cd $(dirname $CUR)/store && pwd)/*.dnn runs/store/ 2>/dev/null
one() {  # NAME IN DEPTH
  D=$H/runs/$1
  # resuming an interrupted level: a chunk that ended with exit 0 and still has its out.state is kept
  [ -s $D/out.state ] && grep -q "^exit 0" $D/cmd.txt 2>/dev/null && return 0
  rm -rf $D; mkdir -p $D
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_WARM_DIR=$H/runs/store
  cd $BNBDIR
  CMD="./bnbdrv -p $H/par/$PAR.par -u $UB -r $2 --allow-settings-change -o $D/out.state -d $3 ../qaplib/$INST.dat"
  { echo "$CMD"; date -u +%FT%TZ; } > $D/cmd.txt
  $CMD > $D/log.txt 2> $D/err.txt; echo "exit $?" >> $D/cmd.txt; date -u +%FT%TZ >> $D/cmd.txt
}
export -f one
for d in $(seq $FROM $UNTIL); do
  [ -e stop ] && { echo "$(date -u +%FT%T) stop file: not starting depth $d" >> progress.txt; break; }
  n=$(grep -c "^node " $CUR); echo "$(date -u +%FT%T) depth $d: $n open nodes" >> progress.txt
  [ "$n" -eq 0 ] && { echo "$(date -u +%FT%T) no open node: the tree is finished" >> progress.txt; break; }
  k=48; [ "$n" -lt "$k" ] && k=$n
  python3 splitstate.py $CUR runs/in_${d}_ $k -contiguous > /dev/null
  seq 0 $((k-1)) | xargs -P 4 -I{} bash -c "one d${d}_c{} $H/runs/in_${d}_{}.state $d"
  bad=$(grep -L "^exit 0" runs/d${d}_c*/cmd.txt | wc -l)
  [ "$bad" -ne 0 ] && { echo "$(date -u +%FT%T) depth $d: $bad processes did not end with exit 0: stopped" >> progress.txt; break; }
  python3 splitstate.py -merge runs/open_$((d+1)).state $(for i in $(seq 0 $((k-1))); do echo runs/d${d}_c$i/out.state; done) > /dev/null
  CUR=runs/open_$((d+1)).state
  python3 openlist.py $CUR > open/open_$((d+1)).txt
  # parents' states that no open node refers to any more
  grep -o "@[^ ]*" $CUR | sort -u | sed 's/^@//' > runs/keep.txt
  for f in runs/store/*.dnn; do grep -q "$(basename $f)" runs/keep.txt || rm -f $f; done
  rm -f runs/in_${d}_*.state runs/d${d}_c*/out.state; [ $d -gt $FROM ] && rm -f runs/open_$d.state
  echo "$(date -u +%FT%T) depth $d done: $(grep -c '^node ' $CUR) open at depth $((d+1)); store $(du -sh runs/store | cut -f1); disk free $(df -h . | tail -1 | awk '{print $4}')" >> progress.txt
  python3 levels.py > logs/levels.out 2>&1; rc=$?
  ( cd /home/user/qap && git add experiments/bound_lab/tai50b_bfs && git commit -q -m "bound_lab/tai50b_bfs: tai50b rehearsal, depth $d complete

Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014vvMYRrFFhnHwcAtHPxwzG" && git push -q origin claude/cool-babbage-ilnogj-tai50b ) >> logs/git.out 2>&1
  [ $rc -ne 0 ] && { echo "$(date -u +%FT%T) levels.py stop condition (exit $rc): see logs/levels.out" >> progress.txt; break; }
done
echo "$(date -u +%FT%T) bfs.sh ended" >> progress.txt
