#!/bin/bash
# chain.sh PREVTAG DEPTH LAST: the on-path chain. From runs/PREVTAG (a node at DEPTH-1 that branched on row r (or column c),
# read from its node trace), bound the child r:p[r] of QAPSOLN/tai60b.sln (onpath.txt) with the incumbent at
# best known + 1 under par/kidW.par, and go on down to depth LAST. Tags dK_R_Vpath.
D=$(cd $(dirname $0); pwd); prev=$1; d=$2; last=$3; OP=${ONPATH:-$D/onpath.txt}; PRE=${PRE:-d}; SUF=${SUF:-path}   # ONPATH/PRE/SUF: another permutation's path (the mirror)
while [ $d -le $last ]; do
  b=$(grep -o " B [0-9]* [rc][0-9]*" $D/runs/$prev/nodes.txt | tail -1 | awk '{print $3}')
  [ -n "$b" ] || { echo "no branching row or column in $prev"; exit 1; }
  if [ "${b:0:1}" = r ]; then r=${b:1}; v=$(tail -1 $OP | tr ' ' '\n' | grep "^$r:" | cut -d: -f2)
  else v=${b:1}; r=$(tail -1 $OP | tr ' ' '\n' | grep ":$v$" | cut -d: -f1); fi   # a column branch: the row holding v
  tag=${PRE}${d}_${r}_${v}${SUF}
  $D/kid.sh $tag kidW $D/runs/$prev/out.state $d $r:$v 608215055
  prev=$tag; d=$((d+1))
done
