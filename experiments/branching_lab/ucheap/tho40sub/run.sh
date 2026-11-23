#!/bin/bash
# ucheap/tho40sub: three depth-3 subtrees of tho40 (from the traced run tho40ab-32: 781, 571 and 466 nodes there),
# each solved to the end as one single-threaded bnbdrv -F run, under tho40's recipe (prod.par =
# experiments/estimates/tho40/tho40.par: rule 9, 2 lines a side) and with the cheap rule (exp_all.par).
# Incumbent 240516 as in the traced run. Each run stops after CAP CPU s (state saved, exit code 3).
#   run.sh        all six, four at a time, longest first
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-cheap}
CAP=${CAP:-4200}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
one() {
  local node=$1 arm=$2 tag; tag=$(echo "$node" | tr ':,' '-_')_$arm
  local t0=$SECONDS
  (cd "$BUILD/latest" && QAP_NODE_TRACE="$D/runs/$tag.nodes.txt" ./bnbdrv -p "$D/$arm.par" -u 240516 -F "$node" -t "$CAP" \
     -o "$D/runs/$tag.state" ../qaplib/tho40t.dat > "$D/runs/$tag.log" 2> "$D/runs/$tag.err"; echo "exit $? wall $((SECONDS - t0)) s" >> "$D/runs/$tag.log")
  echo "$(date -u +%T) $tag: $(grep -E 'Nodes:|Total Time' "$D/runs/$tag.log" | tr -s ' ' | tr '\n' ' ') $(tail -1 "$D/runs/$tag.log")"
}
export -f one; export D BUILD CAP
printf '%s\n' "26:10,8:0,3:32 prod" "26:10,8:0,3:32 exp_all" "26:10,8:0,30:32 prod" "26:10,8:0,36:32 prod" "26:10,8:0,30:32 exp_all" "26:10,8:0,36:32 exp_all" |
  xargs -P 4 -L 1 bash -c 'one $0 $1'
