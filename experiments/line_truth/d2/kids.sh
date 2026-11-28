#!/bin/bash
# d2/kids.sh: a direct check of the pricing at wil50's depth-2 node 21:1,45:0. Each child in kids/list.txt (its
# three fixes; children of c9 and c49, the lines strong branching would take, and of r31, the run's line, as a
# control) is solved to the end from a cold start under the live run's recipe (rule 9 from depth 3 down; no
# firing below depth 1), one process each, four at a time, cheapest by the law first. CAP CPU s a run (exit 3
# with the state saved: the subtree costs at least that). One line a run in kids/done.txt.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); CAP=${CAP:-14400}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
one() {
  local node=$1 tag; tag=$(echo "$node" | tr ':,' '-_')
  grep -q "^$tag " "$D/kids/done.txt" 2>/dev/null && return
  rm -f "$D/kids/$tag.state" "$D/kids/$tag.nodes.txt"
  (cd "$D/../../../latest" && QAP_NODE_TRACE="$D/kids/$tag.nodes.txt" "${BIN:-./bnbdrv}" -p ../infra/phase1/recipes/wil50_d2r10_g150xl2_sb.par \
     -s SYM/wil50t.sym -u 48818 -F "$node" -t "$CAP" -o "$D/kids/$tag.state" ../qaplib/wil50t.dat > "$D/kids/$tag.out.txt" 2>&1
   echo "$tag exit $? $(grep -E "^ *(Nodes:|Total Time)" "$D/kids/$tag.out.txt" | tr -s ' ' | tr '\n' ' ')$(date -u +%FT%TZ)" >> "$D/kids/done.txt")
}
export -f one; export D CAP
awk '{print $1}' "$D/kids/list.txt" | xargs -P 4 -L 1 bash -c 'one $0'
echo "all exit 0 $(date -u +%FT%TZ)" >> "$D/kids/done.txt"
