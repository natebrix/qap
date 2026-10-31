#!/bin/bash
# Elimination rounds on the cloud: elim/round.sh for each node of order.txt,
# 4 at a time. Resumable: a node with a summary_<id>.txt is skipped.
D=$(cd "$(dirname "$0")" && pwd)
one() {
  k=$1; spec=$2
  [ -s "$D/summary_$k.txt" ] && exit 0
  S="$D" "$D/round.sh" "$k" "$spec" > "$D/summary_$k.txt.tmp" 2>&1 && mv "$D/summary_$k.txt.tmp" "$D/summary_$k.txt"
}
export -f one; export D
xargs -P 4 -L 1 bash -c 'one "$@"' _ < "$D/order.txt"
