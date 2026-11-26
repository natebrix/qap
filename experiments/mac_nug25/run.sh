#!/bin/bash
# experiments/mac_nug25/run.sh: nug25's root DNN solve under cognac's solver with this branch's latest/dnn.C
# (the diagnostic for a failed R-step projection, and the stored certificate on an error exit). On an M2
# under Accelerate the solve ended `status error` at iteration 2,662 with bound 1,707 (Nov 26 2026; README.md).
# Builds in a temporary worktree, so no other checkout or running solve is touched. About a minute and a half.
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
# Twice: first as docs/dnn_plus/root_bound.sh ran it (BLAS and LAPACK on one thread: VECLIB_MAXIMUM_THREADS=1,
# where the failure was seen), then with the library's default threading (where it did not fail, Nov 26).
for arm in one-thread default-threads; do
  rm -f "$T/$arm.trace"
  if [ $arm = one-thread ]; then
    VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_DNN_TRACE="$T/$arm.trace" \
      ./bnbdrv -p "$T/plain.par" -u 3746 -d 0 -o "$T/$arm.state" ../qaplib/nug25.dat > "$T/$arm.log" 2>&1 || echo "bnbdrv exit $?"
  else
    QAP_DNN_TRACE="$T/$arm.trace" \
      ./bnbdrv -p "$T/plain.par" -u 3746 -d 0 -o "$T/$arm.state" ../qaplib/nug25.dat > "$T/$arm.log" 2>&1 || echo "bnbdrv exit $?"
  fi
  echo "=== $arm"
  echo "--- diagnostic lines (stderr):"
  grep '^dnn:' "$T/$arm.log" || echo "(none: no projection failed)"
  echo "--- the solve:"
  cat "$T/$arm.trace"
done
echo "--- machine: $(uname -sm)"
echo "(logs: $T/*.log; remove the worktree with: git worktree remove --force $T)"
