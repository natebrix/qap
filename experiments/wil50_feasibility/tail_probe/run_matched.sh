#!/bin/bash
# Run matched tail-probe arms (matched.py) at the same time, half the performance
# cores each, then compare (README.md, "The matched probe"). From anywhere in the repository:
#   experiments/wil50_feasibility/tail_probe/run_matched.sh base full
#   experiments/wil50_feasibility/tail_probe/run_matched.sh baseship fullwide
# Resumable: run the same command again after an interruption.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
[ $# -ge 1 ] || { sed -n 2,6p "$0"; exit 2; }
if [ "$(uname -s)" = Darwin ]; then
  P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
  command -v caffeinate > /dev/null && { caffeinate -ims -w $$ & }
else
  P=$(nproc)
fi
SLOTS=${SLOTS:-$(( P / $# > 1 ? P / $# : 1 ))}
RUNS=${PROBE_RUNS:-$D/runs}     # PROBE_RUNS and PROBE_ARGS: for tests on small instances
cd "$R" || exit 1
[ -x latest/bnbdrv ] && grep -aq QAP_DNN_SHIP latest/bnbdrv || { echo "build first: make solver (needs the shipped-state solver)"; exit 1; }
echo "$(date '+%F %T') arms: $* | $SLOTS slots each | stop with: kill $$ (resume with the same command)"
pids=()
for a in "$@"; do
  mkdir -p "$RUNS"
  python3 "$D/matched.py" "$a" --slots "$SLOTS" --runs "$RUNS" ${PROBE_ARGS:-} > "$RUNS/$a.matched.out" 2>&1 &
  pids+=($!)
done
trap 'for p in "${pids[@]}"; do kill -TERM "$p" 2>/dev/null; done' INT TERM
for p in "${pids[@]}"; do wait "$p"; done
for a in "$@"; do grep -E "^Estimated" "$RUNS/$a.matched.out"; done
python3 "$D/matched.py" --compare --runs "$RUNS"
echo
echo "To share: git add experiments/wil50_feasibility/tail_probe && git commit -m 'matched tail probe: $*' && git push"
