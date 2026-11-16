#!/bin/bash
# One command for the night (bound_lab/laptop/README.md):
#   1. tho40ab_cut.sh   the tho40 subtree with the cut-phase arms, through the coordinator
#   2. nug30.sh         nug30, the record recipe against today's best arm
#   3. tai50b.sh        tai50b whole, if that script is on the branch by then: this step fetches
#                       origin/claude/cool-babbage-ilnogj every 10 minutes for up to WAIT_H hours
#                       (default 6) and starts tai50b.sh from the fetched commit when it is there.
#                       If it never appears, the night ends after step 2.
# From the repository's root, in a worktree of claude/cool-babbage-ilnogj:
#   caffeinate -i experiments/bound_lab/laptop/overnight.sh
# A step that fails does not stop the next one. Everything is logged in overnight.log here.
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../.."
R=$(pwd)
LOG=$H/overnight.log
BRANCH=claude/cool-babbage-ilnogj
WAIT_H=${WAIT_H:-6}
say() { echo "$(date '+%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }

say "overnight: commit $(git rev-parse --short HEAD), $( (sysctl -n hw.ncpu 2>/dev/null || nproc) ) cores"

# 1. tho40ab: the reference arm comes from giveup_ab's results when they are on this machine
ref=""
for c in "$R/experiments/tho40b_analysis/giveup_ab/runs/g200.results.txt" \
         "${COGNAC:-$R/../qap}/experiments/tho40b_analysis/giveup_ab/runs/g200.results.txt"; do
  [ -s "$c" ] && { ref=$c; break; }
done
if [ -n "$ref" ]; then
  say "step 1: tho40ab_cut.sh, arms cut06 comb, reference $ref"
  REF=$ref "$H/tho40ab_cut.sh" >> "$LOG" 2>&1 || say "step 1 failed (exit $?)"
else
  say "step 1: tho40ab_cut.sh, arms g200 cut06 comb (no giveup_ab result found: g200 runs here, about an hour)"
  ARMS="g200 cut06 comb" "$H/tho40ab_cut.sh" >> "$LOG" 2>&1 || say "step 1 failed (exit $?)"
fi

# 2. nug30
say "step 2: nug30.sh"
"$H/nug30.sh" >> "$LOG" 2>&1 || say "step 2 failed (exit $?)"

# 3. tai50b, when its script is on the branch
say "step 3: waiting for experiments/bound_lab/laptop/tai50b.sh on origin/$BRANCH (up to $WAIT_H h)"
tries=$(( WAIT_H * 6 ))
for i in $(seq 1 "$tries"); do
  if git fetch -q origin "$BRANCH" 2>> "$LOG" && git cat-file -e "FETCH_HEAD:experiments/bound_lab/laptop/tai50b.sh" 2> /dev/null; then
    if git checkout -q --detach FETCH_HEAD 2>> "$LOG"; then
      say "step 3: tai50b.sh at commit $(git rev-parse --short HEAD)"
      exec "$H/tai50b.sh" >> "$LOG" 2>&1
    else
      say "step 3: could not check out the fetched commit (local changes?); stopping"
      exit 1
    fi
  fi
  sleep 600
done
say "step 3: tai50b.sh did not appear; done"
