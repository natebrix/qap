#!/bin/bash
# Key off: node traces of the patched build (dom.pair at its default 0) against the
# unpatched build (cognac 0e6c2ecf, bin/bnbdrv.base), on the same searches. Every
# trace field but the last (CPU seconds) must be identical, line for line.
#   bash identity.sh            (one solver process at a time)
D=/home/user/qap/experiments/branching_lab/dom_build
WT=/home/user/wt-dombuild
SCR=/tmp/claude-0/-home-user/92b66c47-7a9f-5ca5-8a65-57e14ac5b075/scratchpad/identity
mkdir -p $SCR $D/logs/identity
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd $WT/latest
one() {  # tag par ub extra... instance
  local tag=$1 par=$2 ub=$3; shift 3
  for b in base patched; do
    local bin=$D/bin/bnbdrv.base; [ $b = patched ] && bin=$WT/latest/bnbdrv
    rm -f $SCR/$tag.$b.trace
    QAP_NODE_TRACE=$SCR/$tag.$b.trace $bin -p $par -u $ub "$@" > $D/logs/identity/$tag.$b.log 2>&1
    echo "$tag $b rc=$? $(grep -m1 'Nodes:' $D/logs/identity/$tag.$b.log)"
  done
  local n=$(wc -l < $SCR/$tag.base.trace)
  if cmp -s <(awk '{$NF=""; print}' $SCR/$tag.base.trace) <(awk '{$NF=""; print}' $SCR/$tag.patched.trace); then
    echo "$tag: IDENTICAL ($n trace lines, all fields but CPU)"
  else
    echo "$tag: DIFFERENT"; diff <(awk '{$NF=""; print}' $SCR/$tag.base.trace) <(awk '{$NF=""; print}' $SCR/$tag.patched.trace) | head -5
  fi
  gzip -c $SCR/$tag.patched.trace > $D/logs/identity/$tag.trace.gz
}
one nug20 $D/par/default.par 2572 ../qaplib/nug20.dat
one tai35b_19-8 $D/par/tai35b_off.par 283315446 -F 19:8 ../qaplib/tai35b.dat
one tho40t_slice $D/par/tho40_off.par 240516 -F 26:10,8:0,16:32,2:8 -s SYM/tho40t.sym ../qaplib/tho40t.dat
