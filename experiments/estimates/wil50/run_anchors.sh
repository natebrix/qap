#!/bin/bash
# wil50: three more complete subtrees of stage-B frontier nodes under the full run's
# settings (fullrun.par), one after another on all performance cores, then the refit.
# From the repository root:
#   experiments/estimates/wil50/run_anchors.sh            # all three, smallest first
#   experiments/estimates/wil50/run_anchors.sh a251       # one of them
# Each is resumable: ^C drains, the same command resumes. Output in
# experiments/estimates/wil50/subtrees/NAME/ (trace.txt, results.txt, summary.txt).
#
# The nodes (stage B's sample, sampling/results/stageB/summary.md; wil50t indices,
# converged shortfall at 4,000 iterations):
#   a220  21:22,31:23,40:15  220.2 short  (the 210-300 class: 30% of the record's total)
#   a251  8:9,21:11,45:40    251.1 short
#   a325  6:49,21:0,45:1     324.6 short  (the tail class: 44%; may be the largest by far)
set -u
R=$(cd "$(dirname "$0")/../../.." && pwd); cd "$R" || exit 1
SET=experiments/wil50_feasibility/fullrun/settings.sh
declare -A SPEC SHORT
SPEC[a220]=21:22,31:23,40:15; SHORT[a220]=220.2
SPEC[a251]=8:9,21:11,45:40;   SHORT[a251]=251.1
SPEC[a325]=6:49,21:0,45:1;    SHORT[a325]=324.6
names=${*:-a220 a251 a325}
for n in $names; do
  [ -n "${SPEC[$n]:-}" ] || { echo "unknown node $n (a220, a251, a325)"; exit 2; }
  EST_RUNS=experiments/estimates/wil50/subtrees experiments/estimates/subtree.sh "$SET" "$n" "${SPEC[$n]}" || exit 1
  [ -f experiments/estimates/wil50/subtrees/$n/trace.txt ] || exit 1   # drained: stop here
done
args=()
for n in a220 a251 a325; do
  [ -f experiments/estimates/wil50/subtrees/$n/results.txt ] && args+=("experiments/estimates/wil50/subtrees/$n=${SHORT[$n]}")
done
python3 experiments/estimates/wil50/refit.py "${args[@]}" | tee experiments/estimates/wil50/refit.txt
echo
echo "To share: git add experiments/estimates/wil50 && git commit -m 'wil50 anchor subtrees' && git push"
