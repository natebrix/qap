#!/bin/bash
# ucheap/sko42sub, second batch: three depth-2 subtrees of sko42 (273, 309 and 310 nodes in sko42-261120b, 7 or 8
# branchings each there, under the depth-1 nodes 36:16, 36:8 and 36:9), as run.sh: one single-threaded bnbdrv -F
# run from the node, incumbent 15814. Arms prod (rule 9), r11x2 (rule 11, two lines a side), r10 (rule 10, every
# line). Each run stops after CAP CPU s (state saved, exit code 3).
#   run_d2.sh        all nine, four at a time; results in runs_d2/, one line a run in runs_d2.log
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
CAP=${CAP:-6600}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$D/runs_d2"
one() {
  local node=$1 arm=$2 tag; tag=$(echo "$node" | tr ':,' '-_')_$arm
  local t0=$SECONDS
  (cd "$BUILD/latest" && QAP_NODE_TRACE="$D/runs_d2/$tag.nodes.txt" ./bnbdrv -p "$D/$arm.par" -s SYM/sko42t.sym -u 15814 -F "$node" -t "$CAP" \
     -o "$D/runs_d2/$tag.state" ../qaplib/sko42t.dat > "$D/runs_d2/$tag.log" 2> "$D/runs_d2/$tag.err"; echo "exit $? wall $((SECONDS - t0)) s" >> "$D/runs_d2/$tag.log")
  echo "$(date -u +%T) $tag: $(grep -E 'Nodes:|Total Time' "$D/runs_d2/$tag.log" | tr -s ' ' | tr '\n' ' ') $(tail -1 "$D/runs_d2/$tag.log")" >> "$D/runs_d2.log"
}
export -f one; export D BUILD CAP
for node in 36:16,23:0 36:8,35:0 36:9,12:0; do for arm in prod r11x2 r10; do echo "$node $arm"; done; done |
  xargs -P 4 -L 1 bash -c 'one $0 $1'
echo done >> "$D/runs_d2.log"
