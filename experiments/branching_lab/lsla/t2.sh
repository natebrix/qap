#!/bin/bash
# lsla T2: a whole tree through the coordinator's local mode (infra/laptop/solve.sh of the BUILD worktree: the
# experiment build plus lsla.patch), as line_signals ran its tai35b / tai40b arms: the params file beside this
# script (line_signals' tai/PROBLEM.par, its region file's path fixed), 2 slots, the incumbent the best known + 1.
#   t2.sh PROBLEM ARM [K=V ...]     e.g. t2.sh tai35b b2 QAP_LS=b:2:0:2
# QAP_LS_PERM (QAPSOLN), QAP_LS_SEARCH (30:3000:2.5:300), QAP_LS_LOG / QAP_LA_LOG (the arm's directory) are set
# here. Outputs in t2/PROBLEM_ARM/: run.log, results.txt, nodes.txt.gz (the node trace, qapctl export),
# ls.log / la.log; the store and work directories are not committed.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-cognac}
P=$1; ARM=$2; shift 2
OUT=$D/t2/${P}_$ARM
mkdir -p "$OUT"
export QAP_LS_PERM=$BUILD/QAPSOLN/$P.sln QAP_LS_SEARCH=${SEARCH:-30:3000:2.5:300} QAP_LS_LOG=$OUT/ls.log QAP_LA_LOG=$OUT/la.log
export QAP_ROOT_CACHE=off
for kv in "$@"; do export "$kv"; done
{ date -u +%FT%TZ; echo "$P $ARM $*"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; md5sum "$BUILD/latest/bnbdrv"; } > "$OUT/cmd.txt"
SLOTS=${SLOTS:-2} OUT=$OUT PARAMS=../qap/experiments/branching_lab/lsla/$P.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
