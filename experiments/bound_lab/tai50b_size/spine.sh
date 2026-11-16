#!/bin/bash
# The path of the best known solution, node by node (validity + how far below the best known value the bound stays):
#   ./spine.sh START.state CORE FROM TO    bounds the on-path node at depths FROM..TO, each warm from its parent with its pool
# (runs/spine_D/: the node at depth D and its children's state). No on-path node may fathom under -u 458821518.
cd "$(dirname "$0")"; IN=$1; CORE=$2
for d in $(seq $3 $4); do
  mkdir -p runs/spine_$d
  python3 splitstate.py $IN runs/spine_$d/in.state -onpath x || exit 1
  ./run_state.sh spine_$d A $PWD/runs/spine_$d/in.state $CORE $d
  IN=runs/spine_$d/out.state
  [ -s $IN ] || exit 0
done
