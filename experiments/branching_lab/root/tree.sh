#!/bin/bash
# root: a grid's whole tree with the root's row forced (Nov 21 2026), to test the root lookahead's ranking.
# Through the coordinator's local mode (infra/laptop/solve.sh of BUILD: cognac 0e6c2ecf + the experiment build
# + ../lsla/lsla.patch + ../tail/tail.patch + root.patch), the nug30 record recipe as it is (rule 9, no cut
# phase, no DNN at the root) with the instance's symmetry file, the incumbent the best known + the step.
#   tree.sh PROBLEM ROW      e.g. tree.sh tho30 25   -> runs/PROBLEM_r25/ (QAP_FORCE_LINE=0:r25)
# The reference is ../grid/runs/PROBLEM_prod (production's own root row; tho30: r5, 14,815 nodes).
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-root}
P=$1; ROW=$2
OUT=$D/runs/${P}_r$ROW
mkdir -p "$OUT"
export QAP_FORCE_LINE=0:r$ROW QAP_ROOT_CACHE=off QAP_LA_LOG=$OUT/la.log
{ date -u +%FT%TZ; echo "$P r$ROW"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; } > "$OUT/cmd.txt"
SLOTS=${SLOTS:-2} OUT=$OUT PARAMS=infra/phase1/recipes/nug30.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
