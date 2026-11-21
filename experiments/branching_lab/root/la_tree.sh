#!/bin/bash
# root: a grid's whole tree with the DNN lookahead at depth 1 and a chosen line score (Nov 21 2026).
#   la_tree.sh PROBLEM ARM [K=V ...]   e.g. la_tree.sh nug30 d1_sum QAP_LA=4:150:1:1
#                                           la_tree.sh nug30 d1_exp QAP_LA=4:150:1:1 QAP_LA_SCORE=exp:0.031
# As tree_tr.sh (the recipe as it is, the node trace on, 2 slots), nothing forced at the root.
# QAP_LA_SCORE: exp:B scores a line by the sum over its children of exp(B * shortfall), the tree they predict
# (B from the production tree's depth-2 nodes: ln(nodes under a node) against its shortfall; nug30 0.031,
# tho30 0.0011); unset, the sum of the shortfalls, as ../grid's runs.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-root}
P=$1; ARM=$2; shift 2
OUT=$D/runs/${P}_$ARM
mkdir -p "$OUT"
export QAP_ROOT_CACHE=off QAP_LA_LOG=$OUT/la.log
for kv in "$@"; do export "$kv"; done
{ date -u +%FT%TZ; echo "$P $ARM $*"; env | grep -E '^QAP_' | sort; git -C "$BUILD" rev-parse --short HEAD; } > "$OUT/cmd.txt"
RC=$(cd "$BUILD" && [ -f infra/phase1/recipes/$P.sh ] && (COORD=; source infra/phase1/recipes/$P.sh > /dev/null; echo "$COORD"))
RC=${RC:---budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300}
case " $RC " in *" --trace "*) ;; *) RC="$RC --trace" ;; esac
echo "COORD=$RC" >> "$OUT/cmd.txt"
COORD="$RC" SLOTS=${SLOTS:-2} OUT=$OUT PARAMS=infra/phase1/recipes/nug30.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
