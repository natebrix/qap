#!/bin/bash
# The matched probe on node 5 in the production configuration (README.md, "Node 5
# in the production configuration"): bestcoret on cores 0-1 and cpcutt on cores 2-3,
# 2 slots each, at the same time, on wil50t.dat with node 5's spec transposed
# (23:49,37:15,6:5 on wil50.dat is 49:23,15:37,5:6), into runs_node5_cpt/
# (PROBE_RUNS to change it). The tasks run latest/bnbdrv: do not rebuild meanwhile.
# Resumable: run it again. The same as
#   PROBE_RUNS=.../runs_node5_cpt PROBE_ARGS="--node 49:23,15:37,5:6 --instance qaplib/wil50t.dat" \
#     SLOTS=2 run_matched.sh bestcoret cpcutt
# with each arm pinned to two cores.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
RUNS=${PROBE_RUNS:-$D/runs_node5_cpt}
ARGS="--node 49:23,15:37,5:6 --instance qaplib/wil50t.dat"
mkdir -p "$RUNS"
cd "$R" || exit 1
[ -x latest/bnbdrv ] && grep -aq "regional cuts only (side 0 or 1)" latest/bnbdrv || {
  echo "build first: make solver (needs the side-1 CP projector)"; exit 1; }
echo "$(date '+%F %T') $(git rev-parse --short HEAD) $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')" >> "$RUNS/hosts.txt"
taskset -c 0,1 python3 "$D/matched.py" bestcoret --slots 2 --runs "$RUNS" $ARGS > "$RUNS/bestcoret.matched.out" 2>&1 &
taskset -c 2,3 python3 "$D/matched.py" cpcutt --slots 2 --runs "$RUNS" $ARGS > "$RUNS/cpcutt.matched.out" 2>&1 &
wait
for a in bestcoret cpcutt; do grep -E "^Estimated" "$RUNS/$a.matched.out"; done
python3 "$D/matched.py" --compare --runs "$RUNS" bestcoret cpcutt
