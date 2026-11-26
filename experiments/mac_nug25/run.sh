#!/bin/bash
# experiments/mac_nug25/run.sh: nug25's root DNN solve under cognac's solver with this branch's latest/dnn.C
# (the diagnostic for a failed R-step projection, and the stored certificate on an error exit). On an M2
# under Accelerate the solve ended `status error` at iteration 2,662 with bound 1,707 (Nov 26 2026; README.md).
# Builds in a temporary worktree, so no other checkout or running solve is touched. About a minute.
# From any checkout of the repository:
#   git fetch origin claude/cool-babbage-ilnogj && git show FETCH_HEAD:experiments/mac_nug25/run.sh | bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
BR=${BRANCH:-claude/cool-babbage-ilnogj}
git -C "$ROOT" fetch -q origin cognac; COG=$(git -C "$ROOT" rev-parse FETCH_HEAD)
git -C "$ROOT" fetch -q origin "$BR"; FIX=$(git -C "$ROOT" rev-parse FETCH_HEAD)
T=${TMPDIR:-/tmp}; T=${T%/}/qap-nug25
git -C "$ROOT" worktree remove --force "$T" 2> /dev/null || rm -rf "$T"
git -C "$ROOT" worktree prune
git -C "$ROOT" worktree add -q --detach "$T" "$COG"
git -C "$ROOT" show "$FIX:latest/dnn.C" > "$T/latest/dnn.C"
git -C "$ROOT" show "$FIX:experiments/mac_nug25/plain.par" > "$T/plain.par"
echo "solver: cognac ${COG:0:8} with latest/dnn.C of $BR ${FIX:0:8}; building in $T"
(cd "$T" && make > build.log 2>&1) || { echo "build failed: $T/build.log"; exit 1; }
cd "$T/latest"
QAP_DNN_TRACE="$T/trace.txt" ./bnbdrv -p "$T/plain.par" -u 3746 -d 0 -o "$T/nug25.state" ../qaplib/nug25.dat > "$T/log.txt" 2>&1 \
  || echo "bnbdrv exit $?"
echo "--- diagnostic lines (stderr):"
grep '^dnn:' "$T/log.txt" || echo "(none: no projection failed)"
echo "--- the solve:"
cat "$T/trace.txt"
echo "--- machine: $(uname -sm), $(grep -m1 -i 'lapack=' "$T/log.txt" || true)"
echo "(log: $T/log.txt; remove the worktree with: git worktree remove --force $T)"
