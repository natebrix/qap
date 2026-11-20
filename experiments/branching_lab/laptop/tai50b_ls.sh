#!/bin/bash
# Overnight on the laptop: tai50b's whole tree with the landscape signal choosing the line at chain-like
# nodes of depths 0-10 (lsla's mode a), everything else tai50b's production recipe (rule 9 at every depth).
#
# The question: on tai50b the signal agrees with rule 9 at depths 0-2 and prefers other lines at depths
# 3-9, with 1-2 predicted survivors against production's 2-6 (../landscape/README.md), unchecked. This is
# the check, on a whole tree. Reference: tai50b-r9, the same recipe without the signal: 6,629 nodes, closed
# at 458,821,517 (cognac experiments/rule9_bottom/README.md). Whole-tree node counts agree between machines
# to about 0.1% (the same README), so nodes are the measure; the two machines' CPU is not comparable.
# Decision rule, fixed before the run: the tree must close at 458,821,517; at most 6,300 nodes (5% fewer)
# and the signal has legs at n = 50; at least 6,960 (5% more) and it is harmful as it stands; between, no
# effect.
#
# From the laptop's qap checkout (any branch; nothing in it is changed):
#
#   git fetch origin cognac claude/cool-babbage-ilnogj
#   git show origin/claude/cool-babbage-ilnogj:experiments/branching_lab/laptop/tai50b_ls.sh > /tmp/tai50b_ls.sh
#   caffeinate -i bash /tmp/tai50b_ls.sh
#
# It makes a worktree beside the checkout (../qap-lsla, cognac 0e6c2ecf), applies the experiment build's
# patch and lsla.patch there, builds, runs a two-minute check on tai25b (the signal must fire and the patched
# solver must run; "CHECK PASSED"), then infra/laptop/solve.sh tai50b in that worktree with the signal on.
# Interrupted: run the same command again (solve.sh resumes from its store). SMOKE_ONLY=1 stops after the
# check. SLOTS as solve.sh's (default: the cores less 2).
set -uo pipefail
R=$(git rev-parse --show-toplevel 2> /dev/null) || { echo "run this inside the qap checkout"; exit 2; }
BASE=0e6c2ecf
SRC=${SRC:-origin/claude/cool-babbage-ilnogj}
W=${W:-$(dirname "$R")/qap-lsla}
P=${PROBLEM:-tai50b}
LS=${LS:-a:2:0:2}
PY=${PYTHON:-python3}
say() { echo "$(date '+%m-%d %H:%M:%S') $*"; }

git -C "$R" cat-file -e "$BASE^{commit}" 2> /dev/null || git -C "$R" fetch -q origin cognac
git -C "$R" cat-file -e "$BASE^{commit}" 2> /dev/null || { echo "commit $BASE (cognac) is not in this checkout: git fetch origin cognac"; exit 2; }
git -C "$R" cat-file -e "$SRC:experiments/branching_lab/lsla/lsla.patch" 2> /dev/null \
  || { echo "no $SRC here: git fetch origin claude/cool-babbage-ilnogj"; exit 2; }

# --- the worktree with the two patches
if [ ! -d "$W" ]; then
  git -C "$R" worktree add -q --detach "$W" "$BASE" || exit 1
  say "worktree $W at $BASE"
fi
cd "$W" || exit 1
[ "$(git rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "$W is not at $BASE"; exit 2; }
if ! grep -q KICKS3 latest/bnb.C; then
  # not patched, or patched with the Nov 19 lsla.patch (before the random-start fix): start latest/ again
  git checkout -q -- latest && git clean -fdq latest
  git apply experiments/line_signals/experiment_build.patch || { echo "experiment_build.patch did not apply"; exit 1; }
  git -C "$R" show "$SRC:experiments/branching_lab/lsla/lsla.patch" > /tmp/lsla.$$.patch
  git apply /tmp/lsla.$$.patch || { echo "lsla.patch did not apply"; exit 1; }
  rm -f /tmp/lsla.$$.patch
  grep -q KICKS3 latest/bnb.C || { echo "$SRC's lsla.patch has no KICKS3: git fetch origin claude/cool-babbage-ilnogj"; exit 1; }
  say "patched: experiment build + lsla (with the random-start fix)"
fi
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# --- build, and the check on tai25b (root and depth 1 with the signal; about 2 minutes)
S=$W/infra/laptop/runs/lsla-check
mkdir -p "$S"
make solver PYTHON="$PY" > "$S/build.log" 2>&1 || { echo "build failed: $S/build.log"; tail -20 "$S/build.log"; exit 1; }
git -C "$R" show "$SRC:experiments/branching_lab/lsla/validate/tai25b_val.par" > "$S/tai25b_val.par"
rm -f "$S/ls.log" "$S/nodes.txt"
( cd latest && QAP_LS=a:2:0:1 QAP_LS_SEARCH=30:0:2.5:300:100 QAP_LS_PERM=../QAPSOLN/tai25b.sln QAP_LS_LOG=$S/ls.log QAP_NODE_TRACE=$S/nodes.txt \
    ./bnbdrv -p "$S/tai25b_val.par" -u 344355647 -d 1 -o "$S/out.state" ../qaplib/tai25b.dat > "$S/log.txt" 2> "$S/err.txt" )
rc=$?
fir=$(grep -c '^ls: ' "$S/ls.log" 2> /dev/null); fir=${fir:-0}
nod=$(grep -c '^node ' "$S/nodes.txt" 2> /dev/null); nod=${nod:-0}
root=$(grep -m1 '^ls: .* depth 0 ' "$S/ls.log" 2> /dev/null | sed 's/.*| chose //')
if [ $rc != 0 ] || [ "$fir" -lt 2 ] || [ "$nod" -lt 20 ] || ! grep -q 'costs 344355646' "$S/err.txt"; then
  echo "CHECK FAILED (exit $rc, $fir signal firings, $nod nodes): $S/err.txt, $S/log.txt"; tail -5 "$S/err.txt"; exit 1
fi
say "CHECK PASSED: tai25b, $fir signal firings, $nod nodes, root line $root (the cloud build: c7)"
[ "${SMOKE_ONLY:-0}" = 1 ] && exit 0

# --- the run
OUT=${OUT:-$W/infra/laptop/runs/$P-ls2}
mkdir -p "$OUT"
export QAP_LS=$LS QAP_LS_PERM=$W/QAPSOLN/$P.sln QAP_LS_SEARCH=${SEARCH:-30:3000:2.5:300:100} QAP_LS_LOG=$OUT/ls.log
export QAP_ROOT_CACHE=off   # the root's line is the signal's: a root kept from another run must not be used
{ date; env | grep -E '^QAP_' | sort; git rev-parse --short HEAD; } >> "$OUT/cmd.txt"
say "$P with QAP_LS=$LS: logs in $OUT (run.log; ls.log, one line per signal firing). The root takes about an hour on one core."
OUT=$OUT infra/laptop/solve.sh "$P"
rc=$?

# --- the summary
echo
echo "=== $P with the landscape signal (QAP_LS=$LS) ==="
[ -s "$OUT/results.txt" ] && grep -E 'Nodes:|Total Time|coordinator:' "$OUT/results.txt"
tr=$OUT/export/nodes.txt.gz
if [ -s "$tr" ]; then
  gzip -dc "$tr" | awk '$1 == "node" {n++; o[$8]++; c += $12; d[$2]++; if ($8 == "B") b[$2]++; if ($2 > m) m = $2}
    END {printf "trace: %d lines (reference tai50b-r9: 6,629 nodes), CPU %.0f s; outcomes", n, c;
         for (k in o) printf " %s %d", k, o[k]; printf "\nnodes / branched by depth:";
         for (i = 0; i <= m && i <= 14; i++) printf " %d:%d/%d", i, d[i], b[i]; printf "\n"}'
fi
if [ -s "$OUT/ls.log" ]; then
  awk '/^ls: /{f++; if ($0 ~ /gate 1/) g++; for (i = 1; i <= NF; i++) if ($i == "cpu_ls") {c += $(i+1); break}}
    END {printf "signal: %d firings, %d past the gate, %.0f CPU s\n", f, g, c}' "$OUT/ls.log"
  echo "lines chosen at depths 0-4:"
  awk -F' [|] ' '/^ls: /{n = split($1, t, " "); dep = -1; fx = "";
      for (i = 1; i <= n; i++) {if (t[i] == "depth") dep = t[i+1]; if (t[i] == "fixes") for (j = i + 1; j <= n; j++) fx = fx " " t[j]}
      if (dep > 4) next; g = ($2 ~ /gate 1$/) ? "passed" : "failed"; ch = "-";
      for (i = 3; i <= NF; i++) if ($i ~ /^chose /) ch = substr($i, 7);
      printf "  depth %d [%s ] gate %s, chose %s\n", dep, fx, g, ch}' "$OUT/ls.log" | head -40
fi
echo "to send back: $OUT/results.txt, $OUT/ls.log, $OUT/export/nodes.txt.gz, $OUT/run.log"
exit $rc
