#!/bin/bash
# chain.sh PREVTAG DEPTH LAST: the on-path chain. From runs/PREVTAG (a node at DEPTH-1 that branched on row r,
# read from its node trace), bound the child r:p[r] of QAPSOLN/tai60b.sln (onpath.txt) with the incumbent at
# best known + 1 under par/kidW.par, and go on down to depth LAST. Tags dK_R_Vpath.
D=$(cd $(dirname $0); pwd); prev=$1; d=$2; last=$3
while [ $d -le $last ]; do
  r=$(grep -o " B [0-9]* r[0-9]*" $D/runs/$prev/nodes.txt | tail -1 | sed 's/.* r//')
  [ -n "$r" ] || { echo "no branching row in $prev"; exit 1; }
  v=$(tail -1 $D/onpath.txt | tr ' ' '\n' | grep "^$r:" | cut -d: -f2)
  tag=d${d}_${r}_${v}path
  $D/kid.sh $tag kidW $D/runs/$prev/out.state $d $r:$v 608215055
  prev=$tag; d=$((d+1))
done
