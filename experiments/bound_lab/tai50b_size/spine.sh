#!/bin/bash
# The path of the best known solution, node by node (validity + how far below the best known value the bound stays):
#   ./spine.sh START.state CORE FROM TO [PERMFILE TAG]  (PERMFILE: follow another permutation, e.g. logs/perm_31_32.txt)
#   ./spine.sh START.state CORE FROM TO    bounds the on-path node at depths FROM..TO, each warm from its parent with its pool
# (runs/spine_D/: the node at depth D and its children's state). No on-path node may fathom under -u 458821518.
cd "$(dirname "$0")"; IN=$1; CORE=$2; PERM=${5:-x}; TAG=${6:-spine}
for d in $(seq $3 $4); do
  mkdir -p runs/${TAG}_$d
  python3 splitstate.py $IN runs/${TAG}_$d/in.state -onpath $PERM || exit 1
  ./run_state.sh ${TAG}_$d A $PWD/runs/${TAG}_$d/in.state $CORE $d
  IN=runs/${TAG}_$d/out.state
  [ -s $IN ] || exit 0
done
