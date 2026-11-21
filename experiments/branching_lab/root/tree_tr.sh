#!/bin/bash
# root: a grid's whole tree with the root's row forced (Nov 21 2026), to test the root lookahead's ranking.
# Through the coordinator's local mode (infra/laptop/solve.sh of BUILD: cognac 0e6c2ecf + the experiment build
# + ../lsla/lsla.patch + ../tail/tail.patch + root.patch), the nug30 record recipe as it is (rule 9, no cut
# phase, no DNN at the root) with the instance's symmetry file, the incumbent the best known + the step.
#   tree_tr.sh PROBLEM ROW (as tree.sh, with the node trace on: runs/PROBLEM_rROW_tr; the trace gives the nodes
#   under each depth-1 child, which expscore.py fits against the lookahead's shortfalls)
#   tree.sh PROBLEM ROW      e.g. tree.sh tho30 25   -> runs/PROBLEM_r25/ (QAP_FORCE_LINE=0:r25)
#   tree.sh PROBLEM prod     the recipe's own root row (nug30: r0), under the same build, slots and incumbent
# tho30's reference is ../grid/runs/tho30_prod (r5, 14,815 nodes, 3 slots).
# Decision rule, fixed before the first tree (Nov 21): over the rows run, the rank correlation between the root
# lookahead's score and the whole tree's nodes; the root lookahead pays on the grid if its best row's tree has
# at most 95% of production's nodes, and its score carries no information there if the correlation is 0.3 or less.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-root}
P=$1; ROW=$2
OUT=$D/runs/${P}_r${ROW}_tr; [ "$ROW" = prod ] && OUT=$D/runs/${P}_prod_tr
mkdir -p "$OUT"
export QAP_ROOT_CACHE=off QAP_LA_LOG=$OUT/la.log
[ "$ROW" = prod ] || export QAP_FORCE_LINE=0:r$ROW      # ROW = prod: the recipe's own root row, nothing forced
{ date -u +%FT%TZ; echo "$P r$ROW"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; } > "$OUT/cmd.txt"
RC=$(cd "$BUILD" && [ -f infra/phase1/recipes/$P.sh ] && (COORD=; source infra/phase1/recipes/$P.sh > /dev/null; echo "$COORD"))
RC=${RC:---budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300}
case " $RC " in *" --trace "*) ;; *) RC="$RC --trace" ;; esac      # the recipe's unit sizing, with the node trace
echo "COORD=$RC" >> "$OUT/cmd.txt"
COORD="$RC" SLOTS=${SLOTS:-2} OUT=$OUT PARAMS=infra/phase1/recipes/nug30.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
