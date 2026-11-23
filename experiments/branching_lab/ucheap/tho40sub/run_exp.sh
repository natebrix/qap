#!/bin/bash
# ucheap/tho40sub: the same three subtrees with the exponential score on rule 9's own lines (exp.par), three at once.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-cheap}
CAP=${CAP:-3000}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for node in 26:10,8:0,3:32 26:10,8:0,30:32 26:10,8:0,36:32; do
  tag=$(echo "$node" | tr ':,' '-_')_exp
  ( t0=$SECONDS; cd "$BUILD/latest" && QAP_NODE_TRACE="$D/runs/$tag.nodes.txt" ./bnbdrv -p "$D/exp.par" -u 240516 -F "$node" -t "$CAP" \
      -o "$D/runs/$tag.state" ../qaplib/tho40t.dat > "$D/runs/$tag.log" 2> "$D/runs/$tag.err"; echo "exit $? wall $((SECONDS - t0)) s" >> "$D/runs/$tag.log"
    echo "$(date -u +%T) $tag: $(grep -E 'Nodes:|Total Time' "$D/runs/$tag.log" | tr -s ' ' | tr '\n' ' ') $(tail -1 "$D/runs/$tag.log")" ) &
done
wait
