#!/bin/bash
# grid: does the DNN lookahead at the top of the tree shrink a grid instance's tree? (Nov 21 2026)
# Whole trees through the coordinator's local mode (infra/laptop/solve.sh of BUILD: cognac 0e6c2ecf + the
# experiment build + lsla.patch), the nug30 record recipe (infra/phase1/recipes/nug30.par: rule 9, no cut phase)
# with the instance's symmetry file, the incumbent the best known + the objective's step.
#   run.sh PROBLEM ARM [K=V ...]     e.g. run.sh tho30 la150 QAP_LA=4:150:0:1
# Outputs in runs/PROBLEM_ARM/: run.log, results.txt, nodes.txt.gz, la.log; stores are not committed.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-lsla-laptop}
P=$1; ARM=$2; shift 2
OUT=$D/runs/${P}_$ARM
mkdir -p "$OUT"
export QAP_LA_LOG=$OUT/la.log QAP_ROOT_CACHE=off
for kv in "$@"; do export "$kv"; done
{ date -u +%FT%TZ; echo "$P $ARM $*"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; } > "$OUT/cmd.txt"
SLOTS=${SLOTS:-3} OUT=$OUT PARAMS=infra/phase1/recipes/nug30.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
