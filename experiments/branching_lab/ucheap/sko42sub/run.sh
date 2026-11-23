#!/bin/bash
# ucheap/sko42sub: three depth-3 subtrees of sko42 (226 or 227 nodes each in the solved run sko42-261120b, under
# three different depth-1 nodes), each solved to the end as one single-threaded bnbdrv -F run from the node (a cold
# first solve there, no inherited pool), incumbent 15814 and the symmetry file as in that run. Arms: prod.par
# (experiments/estimates/sko42/sko42.par: rule 9, one line a side), r11.par (rule 11), r11x2.par (rule 11, two
# lines a side), r10.par (rule 10, every line). Each run stops after CAP CPU s (state saved, exit code 3).
#   run.sh        all twelve, four at a time
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
CAP=${CAP:-4200}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
one() {
  local node=$1 arm=$2 tag; tag=$(echo "$node" | tr ':,' '-_')_$arm
  local t0=$SECONDS
  (cd "$BUILD/latest" && QAP_NODE_TRACE="$D/runs/$tag.nodes.txt" ./bnbdrv -p "$D/$arm.par" -s SYM/sko42t.sym -u 15814 -F "$node" -t "$CAP" \
     -o "$D/runs/$tag.state" ../qaplib/sko42t.dat > "$D/runs/$tag.log" 2> "$D/runs/$tag.err"; echo "exit $? wall $((SECONDS - t0)) s" >> "$D/runs/$tag.log")
  echo "$(date -u +%T) $tag: $(grep -E 'Nodes:|Total Time' "$D/runs/$tag.log" | tr -s ' ' | tr '\n' ' ') $(tail -1 "$D/runs/$tag.log")"
}
export -f one; export D BUILD CAP
for node in 36:1,7:0,5:7 36:16,5:0,22:35 36:3,22:0,24:41; do for arm in prod r11 r11x2 r10; do echo "$node $arm"; done; done |
  xargs -P 4 -L 1 bash -c 'one $0 $1'
