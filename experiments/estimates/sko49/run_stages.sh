#!/bin/bash
# sko49: the wil50 sampling plan (../../wil50_feasibility/sampling/README.md) with
# sko49's current settings (sko49.par, settings.sh), in three steps:
#
#   experiments/estimates/sko49/run_stages.sh a          # stage A: every node at depths 0-2, the depth-3 frontier
#   experiments/estimates/sko49/run_stages.sh b          # stage B: 400 frontier nodes to their give-up, the 80 farthest converged
#   experiments/estimates/sko49/run_stages.sh anchors    # 3 complete subtrees of stage-B nodes, then the estimate
#   experiments/estimates/sko49/run_stages.sh estimate   # the estimate from whatever subtrees have finished
#
# Each step is resumable (run it again) and drains on ^C. On the laptop each takes all
# performance cores; stage A's run directory needs about 10-15 GB for the parents' DNN
# states kept for stage B (tau 0.022: about 430 of an expected 20,000 frontier nodes).
# Output in experiments/estimates/sko49/{stageA,stageB,subtrees}/; stage B's table is
# stageB/summary.md, the estimate estimate.txt.
set -u
STEP=${1:?usage: run_stages.sh a|b|anchors|estimate}
R=$(cd "$(dirname "$0")/../../.." && pwd); cd "$R" || exit 1
H=experiments/estimates/sko49
export SAMPLING_SETTINGS=$R/$H/settings.sh SAMPLING_PARAMS=$R/$H/sko49.par
if [ "$(uname -s)" = Darwin ]; then P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu); else P=$(nproc); fi
CPUS=$(seq -s, 0 $((P - 1)))
S=experiments/wil50_feasibility/sampling
case $STEP in
  a) python3 $S/stage_a.py --run $R/$H/stageA --slots "$P" --cpus "$CPUS" --tau 0.022 ;;
  b) python3 $S/stage_b.py --stage-a $R/$H/stageA --run $R/$H/stageB --n 400 --converge 80 --slots "$P" --cpus "$CPUS" ;;
  anchors|estimate)
    [ -f $H/stageB/summary.md ] || { echo "stage B has not finished ($H/stageB/summary.md)"; exit 1; }
    # three converged misses: the ones nearest the median, the 80th and the 95th percentile of
    # the converged shortfalls (sko49t indices, from stage B's table)
    python3 - $H/stageB/summary.md > $H/anchors.txt <<'PY' || exit 1
import sys
rows = []
for ln in open(sys.argv[1]):
    if not ln.startswith("| ") or ln.startswith("| node") or ln.startswith("|---"): continue
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    if len(c) < 15 or c[5] == "fathomed": continue
    try: rows.append((float(c[12]), c[0]))
    except ValueError: pass
rows.sort()
for name, q in (("q50", 0.50), ("q80", 0.80), ("q95", 0.95)):
    s, spec = rows[min(len(rows) - 1, int(q * len(rows)))]
    print(name, spec, s)
PY
    if [ "$STEP" = anchors ]; then
      while read -r name spec short; do
        EST_RUNS=$H/subtrees experiments/estimates/subtree.sh $H/settings.sh "$name" "$spec" || exit 1
        [ -f $H/subtrees/$name/results.txt ] || exit 1
      done < $H/anchors.txt
    fi
    args=()
    while read -r name spec short; do
      [ -f $H/subtrees/$name/results.txt ] && args+=("$H/subtrees/$name=$short")
    done < $H/anchors.txt
    python3 experiments/estimates/wil50/refit.py --table $H/stageB/summary.md --stage-a $H/stageA --no-record "${args[@]}" | tee $H/estimate.txt ;;
  *) echo "unknown step $STEP (a, b, anchors, estimate)"; exit 2 ;;
esac
echo; echo "To share: git add experiments/estimates/sko49 && git commit -m 'sko49 sampling: step $STEP' && git push"
