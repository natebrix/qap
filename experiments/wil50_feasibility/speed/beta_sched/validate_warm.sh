#!/bin/bash
# Validity of the warm-start key (w=) in the search: tests/brute_force.py with
# warm-started DNN children (-x 0:20:1:200:1, QAP_DNN_WARM_Y=2, as the nug30
# record runs them), for each schedule in SCHEDS, then make dnncheck again
# (switch unset) on the final build. Core 1, one process. Logs: validate/
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
SCHEDS=${SCHEDS:-"two,f=1.6667,t=60,w=0.7 two,f=1.6667,t=5,force=1,w=0.7 w=1.6667"}
mkdir -p "$D/validate"
cd "$R"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
QAP_DNN_WARM_Y=2 taskset -c 1 python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1 > "$D/validate/brute_warm_off.out" 2>&1
echo "brute_force -x 0:20:1:200:1 WARM_Y=2 switch unset: exit $?" >> "$D/validate/summary.txt"
for s in $SCHEDS; do
  tag=${s//[,=]/_}
  QAP_DNN_WARM_Y=2 DNN_BETA_SCHED=$s taskset -c 1 python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1 > "$D/validate/brute_warm_$tag.out" 2>&1
  echo "brute_force -x 0:20:1:200:1 WARM_Y=2 DNN_BETA_SCHED=$s: exit $?" >> "$D/validate/summary.txt"
done
taskset -c 1 make dnncheck > "$D/validate/dnncheck_off_final.out" 2>&1; echo "dnncheck off (final build): exit $?" >> "$D/validate/summary.txt"
