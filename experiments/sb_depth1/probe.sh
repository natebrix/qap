#!/bin/bash
# experiments/sb_depth1/probe.sh [JOBS]: one DNN strong-branching firing at each of sko42's 12 depth-1 nodes
# (README.md): bnbdrv -F 36:k -d 1 under sko42_sb1.par, the node bounded as in sko42-261120b and its line chosen
# by strong branching, then stopped. JOBS firings at a time (default: every core). Resumable: run it again and a
# firing already done is skipped, one under way repeats only the child solves its progress file lacks.
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../.." && pwd)
JOBS=${1:-$(sysctl -n hw.ncpu 2> /dev/null || nproc)}
OUT=$D/runs; mkdir -p "$OUT"
make -C "$R" solver > "$OUT/build.log" 2>&1 || { echo "build failed: $OUT/build.log"; exit 1; }
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 D R OUT
fire() {
  local k=$1 o=$OUT/36_$1
  mkdir -p "$o"
  [ -s "$o/sb.txt" ] && { echo "36:$k done before"; return 0; }
  ( cd "$R/latest" && QAP_SB_TRACE=$o/sb.txt QAP_SB_PROGRESS=$o/progress.txt QAP_DNN_TRACE=$o/dnn.txt \
      ./bnbdrv -p "$D/sko42_sb1.par" -s SYM/sko42t.sym -u 15814 -F "36:$k" -d 1 -o "$o/out.state" \
      "$R/qaplib/sko42t.dat" > "$o/out.txt" 2>&1 )
  echo "36:$k exit $? $(date -u +%H:%M:%S)"
}
export -f fire
echo "$(date -u +%H:%M:%S) 12 firings, $JOBS at a time: $OUT"
# the nine with no symmetry left (rows and columns, about 300 child solves) first, then the three with rows only
printf '%s\n' 9 15 16 8 2 14 1 7 0 10 17 3 | xargs -P "$JOBS" -I{} bash -c 'fire {}'
python3 "$D/summary.py" "$OUT" "$R/infra/phase1/runs/sko42-261120b/nodes.txt.gz" | tee "$OUT/summary.txt"
