#!/bin/bash
# ucheap/tai35b: a whole tree with the experimental lookahead arbiter (branch.la_exp, branch.la_all), through the
# coordinator's local mode, node trace on.
#   run.sh PROBLEM ARM [SLOTS]      e.g. run.sh nug30 exp_all 4     (ARM.par in this directory)
# BUILD (default /home/user/wt-rule10): cognac 7db305d9 (rules 10 and 11), built.
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-rule10}
P=$1; ARM=$2; SLOTS=${3:-4}
OUT=$D/runs/${P}_$ARM
mkdir -p "$OUT" "$BUILD/experiments/branching_lab/ucheap/tai35b"
cp "$D/$ARM.par" "$BUILD/experiments/branching_lab/ucheap/tai35b/$ARM.par"
export QAP_ROOT_CACHE=off
RC=$(cd "$BUILD" && [ -f infra/phase1/recipes/$P.sh ] && (COORD=; source infra/phase1/recipes/$P.sh > /dev/null; echo "$COORD"))
RC=${RC:---budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300}
case " $RC " in *" --trace "*) ;; *) RC="$RC --trace" ;; esac
{ date -u +%FT%TZ; echo "$P $ARM slots $SLOTS"; git -C "$BUILD" rev-parse --short HEAD; echo "COORD=$RC"; } > "$OUT/cmd.txt"
COORD="$RC" SLOTS=$SLOTS OUT=$OUT PARAMS=experiments/branching_lab/ucheap/tai35b/$ARM.par "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
[ -s "$OUT/export/results.txt" ] && cp "$OUT/export/results.txt" "$OUT/results.txt"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
