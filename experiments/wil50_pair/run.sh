#!/bin/bash
# experiments/wil50_pair: rule 10 against rule 9 at depth-2 nodes of wil50, each node's whole subtree under both
# arms, side by side (README.md).
#
#   caffeinate -i experiments/wil50_pair/run.sh [NODE ...]
#
# NODE is a depth-2 node's assignments (-F), e.g. 21:0,45:13; without one, the README's list in its order.
# Each arm runs through experiments/estimates/subtree.sh with half the performance cores (SLOTS=N: another
# number per arm). Ctrl-C drains both arms; the same command goes on where it stopped, and skips what is finished.
# Results: r10/subtrees/NODE/ and r9/subtrees/NODE/ (results.txt, trace.txt, summary.txt), then summary.py.
set -u
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../.." && pwd)
DEFAULT="21:0,45:13 21:3,34:0"
[ $# -gt 0 ] && NODES="$*" || NODES=$DEFAULT
if [ "$(uname -s)" = Darwin ]; then P=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || sysctl -n hw.physicalcpu); else P=$(nproc); fi
SLOTS=${SLOTS:-$(( P / 2 ))}; [ "$SLOTS" -ge 1 ] || SLOTS=1
[ -x "$R/latest/bnbdrv" ] || { echo "no latest/bnbdrv: make -C latest first"; exit 1; }
# Each arm keeps its own local store, whatever the shell says: both arms name the run after the node, so on a
# shared store (EST_BACKEND set, Nov 23) they were one run with two coordinators, closed with nothing merged.
unset EST_BACKEND EST_RUNS
STOP=0; PIDS=""
trap 'STOP=1; echo "stopping: both arms drain (a minute or two); the same command resumes"; [ -z "$PIDS" ] || kill -TERM $PIDS 2> /dev/null' INT TERM
for node in $NODES; do
  [[ $node =~ ^[0-9]+:[0-9]+,[0-9]+:[0-9]+$ ]] || { echo "not a depth-2 node: $node"; exit 2; }
  name=$(echo "$node" | tr ':,' '-_'); PIDS=""
  for arm in r10 r9; do
    A=$D/$arm/subtrees/$name
    if [ -f "$A/init.done" ] && ! grep -q '^local:' "$A/init.done"; then   # made on a shared store: not this arm's own
      mv "$A" "$A.shared-store.$(date +%m%d%H%M%S)"; echo "$node $arm: set aside a run made on a shared store"
    fi
    [ -f "$D/$arm/subtrees/$name/results.txt" ] && { echo "$(date '+%m-%d %H:%M') $node $arm: finished already"; continue; }
    "$R/experiments/estimates/subtree.sh" "$D/$arm/settings.sh" "$name" "$node" "$SLOTS" >> "$D/run.log" 2>&1 &
    PIDS="$PIDS $!"
    echo "$(date '+%m-%d %H:%M') $node $arm: started, $SLOTS slots (log: $D/run.log)"
  done
  for p in $PIDS; do while kill -0 "$p" 2> /dev/null; do wait "$p"; done; done
  [ "$STOP" = 0 ] || break
  for arm in r10 r9; do [ -f "$D/$arm/subtrees/$name/results.txt" ] || { echo "$node $arm did not finish: see $D/run.log and $D/$arm/subtrees/$name/"; exit 1; }; done
  python3 "$D/summary.py" | tee "$D/summary.txt"
done
