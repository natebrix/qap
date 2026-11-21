#!/bin/bash
# Overnight on the laptop: tai50b's whole tree with the rule "tail" (../tail/README.md) in place of rule 9 at
# every node that has DNN marginals; everything else tai50b's production recipe. Nov 21.
#
# The rule: branch on the free line with the fewest marginals above 0.001, ties by rule 9's score. It fixed rule
# 9's four failures on tai60b (0, 0, 0, 2 branching children against 13, 12, 6, 7), picked the best bounded line
# at tai60b's chain nodes at depths 7-9, and cut tai35b's tree from 1,017 to 884 nodes. At tai50b's root it is
# rule 9's own line (every line has 49-50 such entries there), so the root's line should be r31.
#
# The question: does it pay off on a whole tree where rule 9 runs at every depth? Reference: tai50b-r9, the same
# recipe without the rule: 6,629 nodes, closed at 458,821,517 (cognac experiments/rule9_bottom/README.md).
# Whole-tree node counts agree between machines to about 0.1%, so nodes are the measure.
# Decision rule, fixed before the run: the tree must close at 458,821,517; at most 6,300 nodes (5% fewer) and
# the rule helps at n = 50; at least 6,960 (5% more) and it hurts; between, no effect.
# An early reading, about two hours in: the root's line in tail.log (r31 expected) and the branching children at
# depth 1 (10 under r31).
#
# From the laptop's qap checkout (any branch; nothing in it is changed):
#
#   git fetch origin cognac claude/cool-babbage-ilnogj
#   git show origin/claude/cool-babbage-ilnogj:experiments/branching_lab/laptop/tai50b_tail.sh > /tmp/tai50b_tail.sh
#   caffeinate -i bash /tmp/tai50b_tail.sh
#
# It uses the worktree beside the checkout (../qap-lsla, cognac 0e6c2ecf; made if absent), applies the patches
# it lacks (the experiment build, lsla.patch, tail.patch), builds, runs a two-minute check on tai25b ("CHECK
# PASSED"), then infra/laptop/solve.sh tai50b in that worktree, in its own store (infra/laptop/runs/tai50b-tail;
# earlier runs' stores are left alone). Interrupted: run the same command again (solve.sh resumes).
# SMOKE_ONLY=1 stops after the check. TAIL, OUT, SLOTS override the defaults.
set -uo pipefail
R=$(git rev-parse --show-toplevel 2> /dev/null) || { echo "run this inside the qap checkout"; exit 2; }
BASE=0e6c2ecf
SRC=${SRC:-origin/claude/cool-babbage-ilnogj}
W=${W:-$(dirname "$R")/qap-lsla}
P=${PROBLEM:-tai50b}
TAIL=${TAIL:-0.001}
PY=${PYTHON:-python3}
say() { echo "$(date '+%m-%d %H:%M:%S') $*"; }

git -C "$R" cat-file -e "$BASE^{commit}" 2> /dev/null || git -C "$R" fetch -q origin cognac
git -C "$R" cat-file -e "$BASE^{commit}" 2> /dev/null || { echo "commit $BASE (cognac) is not in this checkout: git fetch origin cognac"; exit 2; }
git -C "$R" cat-file -e "$SRC:experiments/branching_lab/tail/tail.patch" 2> /dev/null \
  || { echo "no $SRC here: git fetch origin claude/cool-babbage-ilnogj"; exit 2; }

# --- the worktree with the three patches
if [ ! -d "$W" ]; then
  git -C "$R" worktree add -q --detach "$W" "$BASE" || exit 1
  say "worktree $W at $BASE"
fi
cd "$W" || exit 1
[ "$(git rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "$W is not at $BASE"; exit 2; }
if ! grep -q KICKS3 latest/bnb.C; then
  # not patched, or patched with the Nov 19 lsla.patch: start latest/ again
  git checkout -q -- latest && git clean -fdq latest
  git apply experiments/line_signals/experiment_build.patch || { echo "experiment_build.patch did not apply"; exit 1; }
  git -C "$R" show "$SRC:experiments/branching_lab/lsla/lsla.patch" > /tmp/lsla.$$.patch
  git apply /tmp/lsla.$$.patch || { echo "lsla.patch did not apply"; exit 1; }
  rm -f /tmp/lsla.$$.patch
  say "patched: experiment build + lsla"
fi
if ! grep -q QAP_TAIL latest/bnb.C; then
  git -C "$R" show "$SRC:experiments/branching_lab/tail/tail.patch" > /tmp/tail.$$.patch
  git apply /tmp/tail.$$.patch || { echo "tail.patch did not apply"; exit 1; }
  rm -f /tmp/tail.$$.patch
  say "patched: tail"
fi
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# --- build, and the check on tai25b (root and depth 1 with the rule; about 2 minutes)
S=$W/infra/laptop/runs/tail-check
mkdir -p "$S"
make solver PYTHON="$PY" > "$S/build.log" 2>&1 || { echo "build failed: $S/build.log"; tail -20 "$S/build.log"; exit 1; }
git -C "$R" show "$SRC:experiments/branching_lab/lsla/validate/tai25b_val.par" > "$S/tai25b_val.par"
: > "$S/tail.log"; : > "$S/nodes.txt"
( cd latest && QAP_TAIL=0.001 QAP_TAIL_LOG=$S/tail.log QAP_NODE_TRACE=$S/nodes.txt \
    ./bnbdrv -p "$S/tai25b_val.par" -u 344355647 -d 1 -o "$S/out.state" ../qaplib/tai25b.dat > "$S/log.txt" 2> "$S/err.txt" )
rc=$?
fir=$(grep -c '^tail: ' "$S/tail.log" 2> /dev/null); fir=${fir:-0}
nod=$(grep -c '^node ' "$S/nodes.txt" 2> /dev/null); nod=${nod:-0}
root=$(grep -m1 '^tail: depth 0 ' "$S/tail.log" 2> /dev/null | sed 's/.*| chose //')
got=$(sed 's/.*| chose \([rc][0-9]*\) .*/\1/' "$S/tail.log" | tr '\n' ' ')
if [ $rc != 0 ] || [ "$nod" -lt 20 ] || [ "$got" != "c7 c14 c5 c3 c14 " ]; then
  echo "CHECK FAILED (exit $rc, $fir choices [$got], $nod nodes; the cloud build chose c7 c14 c5 c3 c14): $S/err.txt, $S/log.txt"; tail -5 "$S/err.txt"; exit 1
fi
say "CHECK PASSED: tai25b, $fir choices ($got), $nod nodes, as the cloud build"
[ "${SMOKE_ONLY:-0}" = 1 ] && exit 0

# --- the run
OUT=${OUT:-$W/infra/laptop/runs/$P-tail}
mkdir -p "$OUT"
export QAP_TAIL=$TAIL QAP_TAIL_LOG=$OUT/tail.log
export QAP_ROOT_CACHE=off   # a root kept from another run must not be used
{ date; env | grep -E '^QAP_' | sort; git rev-parse --short HEAD; } >> "$OUT/cmd.txt"
say "$P with QAP_TAIL=$TAIL: logs in $OUT (run.log; tail.log, one line per choice). The root takes about an hour on one core."
OUT=$OUT infra/laptop/solve.sh "$P"
rc=$?

# --- the summary
echo
echo "=== $P with the rule tail (QAP_TAIL=$TAIL) ==="
[ -s "$OUT/results.txt" ] && grep -E 'Nodes:|Total Time|coordinator:' "$OUT/results.txt"
tr=$OUT/export/nodes.txt.gz
if [ -s "$tr" ]; then
  gzip -dc "$tr" | awk '$1 == "node" {n++; o[$8]++; c += $12; d[$2]++; if ($8 == "B") b[$2]++; if ($2 > m) m = $2}
    END {printf "trace: %d lines (reference tai50b-r9: 6,629 nodes), CPU %.0f s; outcomes", n, c;
         for (k in o) printf " %s %d", k, o[k]; printf "\nnodes / branched by depth:";
         for (i = 0; i <= m && i <= 16; i++) printf " %d:%d/%d", i, d[i], b[i]; printf "\n"}'
fi
if [ -s "$OUT/tail.log" ]; then
  echo "choices: $(grep -c '^tail: ' "$OUT/tail.log"); at depths 0-2:"
  awk '/^tail: /{if ($3 <= 2) print "  " $0}' "$OUT/tail.log" | cut -c1-160 | head -30
fi
echo "to send back: $OUT/results.txt, $OUT/tail.log, $OUT/export/nodes.txt.gz, $OUT/run.log"
exit $rc
