#!/bin/bash
# Key on against key off on real subtrees (README, "Savings"): each -F subtree with
# the production recipe (par/*_off.par, *_on.par = the same + dom.pair = 1), both
# arms with the patched build, node traces kept. One solver process at a time.
#   bash slices.sh JOB...      jobs: tai35b6 tho40 tai60b17 tai60b14 tai35b
D=/home/user/qap/experiments/branching_lab/dom_build
WT=/home/user/wt-dombuild
SCR=/tmp/claude-0/-home-user/92b66c47-7a9f-5ca5-8a65-57e14ac5b075/scratchpad/slices
mkdir -p $SCR $D/logs/slices
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd $WT/latest
arm() {  # tag key par ub extra... instance
  local tag=$1 key=$2 par=$3 ub=$4; shift 4
  rm -f $SCR/$tag.$key.trace
  QAP_NODE_TRACE=$SCR/$tag.$key.trace $WT/latest/bnbdrv -p $par -u $ub -o $SCR/$tag.$key.state "$@" > $D/logs/slices/$tag.$key.log 2>&1
  echo "$tag $key rc=$? $(grep -m1 'Nodes:' $D/logs/slices/$tag.$key.log) $(grep 'nodes dropped' $D/logs/slices/$tag.$key.log) $(grep -m1 'Total Time' $D/logs/slices/$tag.$key.log)"
  gzip -c $SCR/$tag.$key.trace > $D/logs/slices/$tag.$key.trace.gz
}
for job in "$@"; do case $job in
  tai35b6)  # tai35b: depth-6 node with 14 of its 58 trace nodes under dominated nodes
    F=19:26,18:12,23:5,26:6,27:34,30:14
    arm tai35b6 off $D/par/tai35b_off.par 283315446 -F $F ../qaplib/tai35b.dat
    arm tai35b6 on  $D/par/tai35b_on.par  283315446 -F $F ../qaplib/tai35b.dat ;;
  tho40)    # tho40t: the identity slice (175 trace nodes, 14 under dominated nodes); off = identity.sh's run
    arm tho40 on $D/par/tho40_on.par 240516 -F 26:10,8:0,16:32,2:8 -s SYM/tho40t.sym ../qaplib/tho40t.dat ;;
  tai60b17) # tai60b: path chain depth-17 node, 44 trace nodes, 43 under dominated nodes
    F=44:45,43:58,1:33,6:24,7:19,0:32,10:59,14:54,19:2,8:29,4:23,2:49,5:57,9:55,12:0,3:36,11:39
    arm tai60b17 off $D/par/tai60b_off.par 608215055 -F $F ../qaplib/tai60b.dat
    arm tai60b17 on  $D/par/tai60b_on.par  608215055 -F $F ../qaplib/tai60b.dat ;;
  tai60b14) # tai60b: path chain depth-14 node, 179 trace nodes (0.86 CPU h), 169 under dominated nodes
    F=44:45,43:58,1:33,6:24,7:19,0:32,10:39,14:54,19:2,8:29,4:23,2:49,5:55,9:57
    arm tai60b14 on  $D/par/tai60b_on.par  608215055 -F $F ../qaplib/tai60b.dat
    arm tai60b14 off $D/par/tai60b_off.par 608215055 -F $F ../qaplib/tai60b.dat ;;
  tai35b)   # tai35b whole, key on (the key-off tree: line_signals/runs/tai35b_base)
    arm tai35b on $D/par/tai35b_on.par 283315446 ../qaplib/tai35b.dat ;;
esac; done
