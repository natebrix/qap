#!/bin/bash
# The tho40 method (../tho40_estimate/README.md) for any instance with a settings file
# here: (1) the production search, traced, for HOURS of wall; (2) the depth-2 sweep,
# every depth-1 node bounded with its children (bnbdrv -F <node> -d 2), so nothing above
# depth 2 is extrapolated; (3) ../tho40_estimate/estimate.py: log10(subtree size) fitted
# to the gap per depth on the finished subtrees, applied to the open nodes.
#
#   experiments/estimates/trace_estimate.sh INSTANCE [HOURS] [JOBS]
#   experiments/estimates/trace_estimate.sh tho40 5        # the Nov 8 protocol: 5 h
#   experiments/estimates/trace_estimate.sh sko42 6
#
# INSTANCE is a directory here with settings.sh (tho40, sko42). JOBS defaults to the
# performance cores. Each step is skipped when its output exists, so an interrupted run
# resumes with the same command (the search from its checkpoint). Output in
# experiments/estimates/INSTANCE/trace/: search.trace, depth2.trace, estimate.txt; the
# checkpoint (large) in search.ckpt*, ignored by git.
set -u
[ $# -ge 1 ] || { sed -n '2,/^set -u$/p' "$0" | sed '$d'; exit 2; }
INST=$1; HOURS=${2:-5}
R=$(cd "$(dirname "$0")/../.." && pwd); cd "$R" || exit 1
. "experiments/estimates/$INST/settings.sh" || exit 1
fullrun_check || exit 1
if [ "$(uname -s)" = Darwin ]; then
  P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
  command -v caffeinate > /dev/null && { caffeinate -ims -w $$ & }
else P=$(nproc); fi
JOBS=${3:-$P}
O=${EST_OUT:-experiments/estimates/$INST/trace}; mkdir -p "$O/depth2"   # EST_OUT: another output directory
OUT=$(cd "$O" && pwd)
N=$(head -1 "$FULLRUN_INSTANCE" | tr -d ' \r')
SYMF=(); [ -n "$FULLRUN_SYM" ] && SYMF=(-s "$FULLRUN_SYM")
log() { echo "$(date '+%F %T') $INST: $*" | tee -a "$OUT/run.log"; }

# (1) the search, traced, with a wall limit; a resume continues from the checkpoint
if [ ! -f "$OUT/search.done" ]; then
  RES=(); [ -f "$OUT/search.ckpt" ] && RES=(--resume "$OUT/search.ckpt")
  log "search: $JOBS jobs, ${HOURS} h of wall${RES:+ (resumed)}"
  QAP_DNN_TRACE=$OUT/search.trace python3 bench/parallel.py "$FULLRUN_INSTANCE" --jobs "$JOBS" --ub "$FULLRUN_UB" \
    --params "$FULLRUN_PAR" ${FULLRUN_SYM:+--extra "-s $FULLRUN_SYM"} --checkpoint "$OUT/search.ckpt" --checkpoint-every 3600 \
    --wall-limit "$(python3 -c "print(int($HOURS*3600))")" "${RES[@]}" > "$OUT/search.log" 2>&1
  rc=$?
  case $rc in
    0) log "search finished: the tree is solved"; touch "$OUT/search.done" "$OUT/search.solved" ;;
    3) log "search stopped at the wall limit"; touch "$OUT/search.done" ;;
    *) log "search exited $rc (see $O/search.log); run again to resume"; exit 1 ;;
  esac
fi

# (2) the depth-2 sweep: the root's children from a depth-0 pass, then one depth-limited run each
if [ ! -f "$OUT/depth2.done" ]; then
  ( cd latest && ./bnbdrv -p "$FULLRUN_PAR" -u "$FULLRUN_UB" "${SYMF[@]}" -d 0 -o "$OUT/depth2/root.state" "$FULLRUN_INSTANCE" > "$OUT/depth2/root.log" 2>&1 ) \
    || { log "the depth-0 pass failed (see $O/depth2/root.log)"; exit 1; }
  python3 - "$OUT/depth2/root.state" "$N" > "$OUT/depth2/nodes.txt" <<'PY' || exit 1
import sys
sys.path.insert(0, "bench")
from parallel import State
st = State(sys.argv[1]); n = int(sys.argv[2])
for text, _ in st.nodes:
    p = [int(x) for x in text.split()[5:5 + n]]
    print(",".join(f"{i}:{j}" for i, j in enumerate(p) if j >= 0))
PY
  # EST_SWEEP_LIMIT=K: the first K depth-1 nodes only (a smoke test; the sweep stays unfinished)
  [ -n "${EST_SWEEP_LIMIT:-}" ] && { head -n "$EST_SWEEP_LIMIT" "$OUT/depth2/nodes.txt" > "$OUT/depth2/nodes.run"; } || cp "$OUT/depth2/nodes.txt" "$OUT/depth2/nodes.run"
  log "depth-2 sweep: $(wc -l < "$OUT/depth2/nodes.txt") depth-1 nodes, $(wc -l < "$OUT/depth2/nodes.run") to run, $JOBS at a time"
  export FULLRUN_PAR FULLRUN_UB FULLRUN_INSTANCE OUT
  export SYMARG="${FULLRUN_SYM:+-s $FULLRUN_SYM}"
  xargs -P "$JOBS" -I{} sh -c 'f={}; t=$(echo "$f" | tr ":," "_-"); [ -f "$OUT/depth2/d2_$t.done" ] && exit 0;
      cd latest && QAP_DNN_TRACE="$OUT/depth2/d2_$t.trace" ./bnbdrv -p "$FULLRUN_PAR" -u "$FULLRUN_UB" $SYMARG -F "$f" -d 2 "$FULLRUN_INSTANCE" \
      > "$OUT/depth2/d2_$t.log" 2>&1 && touch "$OUT/depth2/d2_$t.done"' < "$OUT/depth2/nodes.run"
  n1=$(wc -l < "$OUT/depth2/nodes.txt"); nd=$(ls "$OUT"/depth2/d2_*.done 2>/dev/null | wc -l)
  [ "$nd" -eq "$n1" ] || { log "depth-2 sweep: $nd of $n1 done; run again to finish"; exit 1; }
  cat "$OUT"/depth2/d2_*.trace > "$OUT/depth2.trace"; touch "$OUT/depth2.done"
  log "depth-2 sweep done: $(grep -c ' depth 2 ' "$OUT/depth2.trace") depth-2 nodes"
fi

# (3) the estimate, self-calibrated and calibrated on the solved tho30 and nug30 trees
N1=$(wc -l < "$OUT/depth2/nodes.txt")
T=$OUT/search.trace+$OUT/depth2.trace
{
  echo "# $INST: trace estimate ($(date '+%F'); search $(grep -c '^pid ' "$OUT/search.trace") traced nodes, depth-2 sweep over $N1 depth-1 nodes)"
  for cal in self tho30 nug30; do
    echo; echo "## calibrated on: $cal"
    case $cal in
      self)  EXHAUSTIVE=1 python3 -I experiments/tho40_estimate/estimate.py predict "$T" "$N" "$N1" "$T" "$N" "$N1" ;;
      tho30) EXHAUSTIVE=1 python3 -I experiments/tho40_estimate/estimate.py predict <(gzip -dc experiments/tho40_estimate/traces/tho30.trace.gz) 30 10 "$T" "$N" "$N1" ;;
      nug30) EXHAUSTIVE=1 python3 -I experiments/tho40_estimate/estimate.py predict <(gzip -dc experiments/tho40_estimate/traces/nug30.trace.gz) 30 9 "$T" "$N" "$N1" ;;
    esac
  done
} > "$OUT/estimate.txt" 2>&1
log "estimate written to $O/estimate.txt"
tail -n 40 "$OUT/estimate.txt"
echo; echo "To share: git add experiments/estimates/$INST && git commit -m '$INST trace estimate' && git push"
