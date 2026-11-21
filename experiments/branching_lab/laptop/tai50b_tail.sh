#!/bin/bash
# Overnight on the laptop: tai50b's whole tree with the gated rule "tail" (../tail/README.md) in place of rule 9
# at the nodes where the rule's count is small; everything else tai50b's production recipe. Nov 21, second version.
#
# The rule: branch on the free line with the fewest marginals above 0.001, ties by rule 9's score; with the gate
# (QAP_TAIL_GATE=0.5) only where that count is at most half the node's free dimension, rule 9 elsewhere.
# What is known: one level, it fixed rule 9's four failures on tai60b (0, 0, 0, 2 branching children against
# 13, 12, 6, 7) and picked the best bounded line at tai60b's chain nodes at depths 7-9. On whole trees the
# ungated rule helped on tai35b (884 nodes against 1,017) and hurt on tai40b (a flat root: at least 843 against
# 731), and the gated rule had no effect on either (1,045; 727): those trees have no tight node going wrong
# under rule 9 for it to repair. tai50b's tree is the first here that may. Its root is flat (49-50 entries on
# every line), so the root is gated and its line is rule 9's r31.
#
# The question: does it pay off on a whole tree where rule 9 runs at every depth? Reference: tai50b-r9, the same
# recipe without the rule: 6,629 nodes, closed at 458,821,517 (cognac experiments/rule9_bottom/README.md).
# Whole-tree node counts agree between machines to about 0.1%, so nodes are the measure.
# Decision rule, fixed before the run: the tree must close at 458,821,517; at most 6,300 nodes (5% fewer) and
# the rule helps at n = 50; at least 6,960 (5% more) and it hurts; between, no effect.
# An early reading, about two hours in: the root "gated" in tail.log, r31 in the trace, 10 branching children at
# depth 1.
#
# From the laptop's qap checkout (any branch; nothing in it is changed):
#
#   git fetch origin cognac claude/cool-babbage-ilnogj
#   git show origin/claude/cool-babbage-ilnogj:experiments/branching_lab/laptop/tai50b_tail.sh > /tmp/tai50b_tail.sh
#   caffeinate -i bash /tmp/tai50b_tail.sh
#
# It uses the worktree beside the checkout (../qap-lsla, cognac 0e6c2ecf; made if absent), applies the patches
# it lacks (the experiment build, lsla.patch, tail.patch), builds, runs a two-minute check on tai25b ("CHECK
# PASSED"), then infra/laptop/solve.sh tai50b in that worktree, in its own store (infra/laptop/runs/tai50b-tailg;
# earlier runs' stores are left alone). Interrupted: run the same command again (solve.sh resumes).
# SMOKE_ONLY=1 stops after the check. TAIL, GATE, OUT, SLOTS override the defaults.
set -uo pipefail
R=$(git rev-parse --show-toplevel 2> /dev/null) || { echo "run this inside the qap checkout"; exit 2; }
BASE=0e6c2ecf
SRC=${SRC:-origin/claude/cool-babbage-ilnogj}
W=${W:-$(dirname "$R")/qap-lsla}
P=${PROBLEM:-tai50b}
TAIL=${TAIL:-0.001}
GATE=${GATE:-0.5}
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
if ! grep -q KICKS3 latest/bnb.C || { grep -q QAP_TAIL latest/bnb.C && ! grep -q QAP_TAIL_GATE latest/bnb.C; }; then
  # not patched, or patched with the Nov 19 lsla.patch or the tail.patch before the gate: start latest/ again
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
  grep -q QAP_TAIL_GATE latest/bnb.C || { echo "$SRC's tail.patch has no gate: git fetch origin claude/cool-babbage-ilnogj"; exit 1; }
  say "patched: tail (with the gate)"
fi
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# --- build, and the check on tai25b (root and depth 1 with the rule; about 2 minutes)
S=$W/infra/laptop/runs/tail-check
mkdir -p "$S"
make solver PYTHON="$PY" > "$S/build.log" 2>&1 || { echo "build failed: $S/build.log"; tail -20 "$S/build.log"; exit 1; }
git -C "$R" show "$SRC:experiments/branching_lab/lsla/validate/tai25b_val.par" > "$S/tai25b_val.par"
: > "$S/tail.log"; : > "$S/nodes.txt"
( cd latest && QAP_TAIL=0.001 QAP_TAIL_GATE=0.5 QAP_TAIL_LOG=$S/tail.log QAP_NODE_TRACE=$S/nodes.txt \
    ./bnbdrv -p "$S/tai25b_val.par" -u 344355647 -d 1 -o "$S/out.state" ../qaplib/tai25b.dat > "$S/log.txt" 2> "$S/err.txt" )
rc=$?
fir=$(grep -c '^tail: ' "$S/tail.log" 2> /dev/null); fir=${fir:-0}
nod=$(grep -c '^node ' "$S/nodes.txt" 2> /dev/null); nod=${nod:-0}
want="gated:c7 gated:c21 gated:c14 chose:c7 "
got=$(sed 's/.*| \([a-z]*\) \([rc][0-9]*\) .*/\1:\2/' "$S/tail.log" | tr '\n' ' ')
if [ $rc != 0 ] || [ "$nod" -lt 20 ] || [ "$got" != "$want" ]; then
  echo "CHECK FAILED (exit $rc, $fir lines [$got], $nod nodes; the cloud build: $want): $S/err.txt, $S/log.txt"; tail -5 "$S/err.txt"; exit 1
fi
say "CHECK PASSED: tai25b, $fir choices ($got), $nod nodes, as the cloud build"
[ "${SMOKE_ONLY:-0}" = 1 ] && exit 0

# --- the run
OUT=${OUT:-$W/infra/laptop/runs/$P-tailg}
mkdir -p "$OUT"
export QAP_TAIL=$TAIL QAP_TAIL_GATE=$GATE QAP_TAIL_LOG=$OUT/tail.log
export QAP_ROOT_CACHE=off   # a root kept from another run must not be used
{ date; env | grep -E '^QAP_' | sort; git rev-parse --short HEAD; } >> "$OUT/cmd.txt"
say "$P with QAP_TAIL=$TAIL QAP_TAIL_GATE=$GATE: logs in $OUT (run.log; tail.log, one line per choice). The root takes about an hour on one core."
OUT=$OUT infra/laptop/solve.sh "$P"
rc=$?

# --- the summary
echo
echo "=== $P with the gated rule tail (QAP_TAIL=$TAIL QAP_TAIL_GATE=$GATE) ==="
[ -s "$OUT/results.txt" ] && grep -E 'Nodes:|Total Time|coordinator:' "$OUT/results.txt"
tr=$OUT/export/nodes.txt.gz
if [ -s "$tr" ]; then
  gzip -dc "$tr" | awk '$1 == "node" {n++; o[$8]++; c += $12; d[$2]++; if ($8 == "B") b[$2]++; if ($2 > m) m = $2}
    END {printf "trace: %d lines (reference tai50b-r9: 6,629 nodes), CPU %.0f s; outcomes", n, c;
         for (k in o) printf " %s %d", k, o[k]; printf "\nnodes / branched by depth:";
         for (i = 0; i <= m && i <= 16; i++) printf " %d:%d/%d", i, d[i], b[i]; printf "\n"}'
fi
if [ -s "$OUT/tail.log" ]; then
  echo "nodes with marginals: $(grep -c '^tail: ' "$OUT/tail.log"), the rule's line at $(grep -c '| chose ' "$OUT/tail.log"), gated (rule 9) at $(grep -c '| gated ' "$OUT/tail.log"); at depths 0-2:"
  awk '/^tail: /{if ($3 <= 2) print "  " $0}' "$OUT/tail.log" | cut -c1-160 | head -30
fi
echo "to send back: $OUT/results.txt, $OUT/tail.log, $OUT/export/nodes.txt.gz, $OUT/run.log"
exit $rc
