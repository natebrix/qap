#!/bin/bash
# The rule "tail" on a whole tree, through the coordinator's local mode (as ../lsla/t2.sh): the params file of
# ../lsla (line_signals' tai/PROBLEM.par), 2 slots, the incumbent the best known + 1, the root not from a cache.
#   run.sh PROBLEM ARM [K=V ...]     e.g. run.sh tai35b tail QAP_TAIL=0.001
# BUILD: a worktree of cognac 0e6c2ecf with experiments/line_signals/experiment_build.patch, ../lsla/lsla.patch
# and tail.patch applied, built with `make solver`. QAP_TAIL=THR[:DMIN:DMAX]: at a node with DNN marginals
# (depths DMIN..DMAX; default all) branch on the free line with the fewest marginals above THR, ties by the larger
# sum of squares, then rows before columns. In place of the params' rule (rule 9 through depth 6, rule 2 below).
# Outputs in runs/PROBLEM_ARM/: run.log, results.txt, nodes.txt.gz, tail.log (one line per choice).
#
# Decision rule, fixed before the first run (Nov 21). Reference: ../lsla/t2/PROBLEM_prod, the same params without
# the rule: tai35b 1,017 nodes, 6,935 CPU s; tai40b 731 nodes. The tree must close at the best known value.
# At most 95% of the reference's nodes: the rule helps on a whole tree. At least 105%: it hurts. Between: no effect.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-tail}
P=$1; ARM=$2; shift 2
OUT=$D/runs/${P}_$ARM
mkdir -p "$OUT"
grep -q QAP_TAIL "$BUILD/latest/bnb.C" || { echo "$BUILD/latest has no tail.patch"; exit 2; }
export QAP_TAIL_LOG=$OUT/tail.log QAP_ROOT_CACHE=off
for kv in "$@"; do export "$kv"; done
{ date -u +%FT%TZ; echo "$P $ARM $*"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; } > "$OUT/cmd.txt"
SLOTS=${SLOTS:-2} OUT=$OUT PARAMS=../qap/experiments/branching_lab/lsla/$P.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
