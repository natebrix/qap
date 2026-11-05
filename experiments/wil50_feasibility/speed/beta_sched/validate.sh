#!/bin/bash
# Validity checks on core 1, one process at a time. Logs: validate/*.out
#  1. make dnncheck (tests/dnn_test.py, all parts), switch off
#  2. tests/dnn_test.py parts A, C, D, E with the switch on (part B compares
#     dnntest -P with the fixed-beta prototype dnn.py, so it does not apply),
#     for each schedule in SCHEDS
#  3. tests/brute_force.py --solver-args "-x 0:20" (the DNN bound in the
#     search) with the switch on, for each schedule in SCHEDS
# SCHEDS: the DNN_BETA_SCHED values tested. t=5 with force=1 makes every
# node of more than 5 iterations change beta.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
SCHEDS=${SCHEDS:-"two,f=1.6667,t=60 two,f=1.6667,t=5,force=1 two,f=0.6667,t=5,force=1 rb,c=8,every=10"}
mkdir -p "$D/validate"
cd "$R"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [ -z "$SKIP_DEFAULT" ]; then
  taskset -c 1 make dnncheck > "$D/validate/dnncheck_off.out" 2>&1; echo "dnncheck off: exit $?" >> "$D/validate/summary.txt"
fi
for s in $SCHEDS; do
  tag=${s//[,=]/_}
  DNN_BETA_SCHED=$s taskset -c 1 python3 tests/dnn_test.py --parts ACDE > "$D/validate/dnntest_$tag.out" 2>&1
  echo "dnn_test ACDE DNN_BETA_SCHED=$s: exit $?" >> "$D/validate/summary.txt"
  DNN_BETA_SCHED=$s taskset -c 1 python3 tests/brute_force.py --solver-args "-x 0:20" --jobs 1 > "$D/validate/brute_$tag.out" 2>&1
  echo "brute_force -x 0:20 DNN_BETA_SCHED=$s: exit $?" >> "$D/validate/summary.txt"
done
