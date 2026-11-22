#!/bin/bash
# port: a whole tree with cognac's DNN strong branching (sb.*), through the coordinator's local mode.
#   run.sh PROBLEM ARM PAR [SLOTS]      e.g. run.sh nug30 sb nug30_sb.par 4
# PAR is a params file in this directory; BUILD (default /home/user/wt-port) is a checkout of cognac
# with the port applied and built. The node trace is on; QAP_SB_TRACE and QAP_SB_PROGRESS go to the
# run's directory (sb.txt, sb_progress.txt).
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-port}
P=$1; ARM=$2; PAR=$3; SLOTS=${4:-4}
OUT=$D/runs/${P}_$ARM
mkdir -p "$OUT" "$BUILD/experiments/branching_lab/port"
cp "$D/$PAR" "$BUILD/experiments/branching_lab/port/$PAR"
export QAP_ROOT_CACHE=off QAP_SB_TRACE=$OUT/sb.txt QAP_SB_PROGRESS=$OUT/sb_progress.txt
RC=$(cd "$BUILD" && [ -f infra/phase1/recipes/$P.sh ] && (COORD=; source infra/phase1/recipes/$P.sh > /dev/null; echo "$COORD"))
RC=${RC:---budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300}
case " $RC " in *" --trace "*) ;; *) RC="$RC --trace" ;; esac
{ date -u +%FT%TZ; echo "$P $ARM $PAR slots $SLOTS"; git -C "$BUILD" rev-parse --short HEAD; echo "COORD=$RC"; } > "$OUT/cmd.txt"
COORD="$RC" SLOTS=$SLOTS OUT=$OUT PARAMS=experiments/branching_lab/port/$PAR "$BUILD/infra/laptop/solve.sh" "$P"
rc=$?
[ -s "$OUT/export/nodes.txt.gz" ] && cp "$OUT/export/nodes.txt.gz" "$OUT/nodes.txt.gz"
echo "exit $rc $(date -u +%FT%TZ)" >> "$OUT/cmd.txt"
exit $rc
