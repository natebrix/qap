#!/bin/bash
# ucheap/sko42sub, third batch: large depth-2 subtrees of sko42 (the 27 largest of the 269 hold 61% of
# sko42-261120b's nodes), as run.sh: one single-threaded bnbdrv -F run from the node, incumbent 15814.
#   run_big.sh NODE ARM...      e.g. run_big.sh 36:14,7:0 prod r10     (the arms side by side; results in runs_big/)
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
CAP=${CAP:-16000}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$D/runs_big"
node=$1; shift
for arm in "$@"; do
  tag=$(echo "$node" | tr ':,' '-_')_$arm
  ( t0=$SECONDS; cd "$BUILD/latest" && QAP_NODE_TRACE="$D/runs_big/$tag.nodes.txt" ./bnbdrv -p "$D/$arm.par" -s SYM/sko42t.sym -u 15814 -F "$node" -t "$CAP" \
      -o "$D/runs_big/$tag.state" ../qaplib/sko42t.dat > "$D/runs_big/$tag.log" 2> "$D/runs_big/$tag.err"; echo "exit $? wall $((SECONDS - t0)) s" >> "$D/runs_big/$tag.log"
    echo "$(date -u +%T) $tag: $(grep -E 'Nodes:|Total Time' "$D/runs_big/$tag.log" | tr -s ' ' | tr '\n' ' ') $(tail -1 "$D/runs_big/$tag.log")" >> "$D/runs_big.log" ) &
done
wait
echo "done $node" >> "$D/runs_big.log"
