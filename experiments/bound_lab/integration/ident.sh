#!/bin/bash
# ident.sh NAME BUILD ARGS... : one traced bnbdrv run of one build, for the switch-off identity.
#   BUILD base = the unpatched tree (worktree /home/user/wt-base, this branch's HEAD before the merge,
#   latest/ identical to cognac e1596e22); new = the merged tree (/home/user/qap).
# Output ident/NAME_BUILD/{cmd.txt,log.txt,err.txt,nodes.txt,dnn.txt}. Paths in ARGS are relative to latest/.
# The binary is copied to bin_local/ (not committed) first, so a rebuild during a long run cannot touch it;
# its sha256 is recorded in cmd.txt.
set -u
H=$(cd "$(dirname "$0")" && pwd); NAME=$1; BUILD=$2; shift 2
case $BUILD in base) R=/home/user/wt-base;; new) R=/home/user/qap;; *) echo "BUILD base|new"; exit 2;; esac
D=$H/ident/${NAME}_$BUILD; mkdir -p "$D" "$H/bin_local"
BIN=$H/bin_local/bnbdrv_${BUILD}_$(sha256sum "$R/latest/bnbdrv" | cut -c1-12)
[ -x "$BIN" ] || cp "$R/latest/bnbdrv" "$BIN"
rm -f "$D/nodes.txt" "$D/dnn.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt
cd "$R/latest"
{ echo "cwd $R/latest"; echo "$(basename "$BIN") $*"; env | grep -E '^(QAP_|DNN_)' | sed "s|$H/||" | sort; date -u +%FT%TZ; } > "$D/cmd.txt"
"$BIN" "$@" > "$D/log.txt" 2> "$D/err.txt"
echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt"
