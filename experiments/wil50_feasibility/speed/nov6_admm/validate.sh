#!/bin/bash
# Validation of the DNN_ADMM_* switches and of the default path, on core CPU (1):
#   1. make test, make dnncheck, make shipcheck, make continuationcheck and
#      python3 bench/bench.py check with every switch unset (the default path);
#   2. tests/dnn_test.py --parts ACDE (bound validity against brute force on
#      2,592 node bounds, fathoming on the nug20 study nodes) with
#      DNN_ADMM_CROSS=2, then with the two warm switches;
#   3. tests/brute_force.py with -x 0:20:1:200:1 and QAP_DNN_WARM_Y=2 (warm
#      children) under the three switches together;
#   4. with BASE set (the dnntest built before the change): identity_check.sh,
#      the default path against that binary.
# Logs: validate/*.log; the last line of each is its verdict.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CPU=${CPU:-1}
mkdir -p "$D/validate"
cd "$R" || exit 1
unset DNN_ADMM_CROSS DNN_ADMM_WARM_BETA DNN_ADMM_WARM_DAMP DNN_BETA_SCHED
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run() { local name=$1; shift; echo "== $name: $*"; taskset -c "$CPU" "$@" > "$D/validate/$name.log" 2>&1; echo "   exit $? : $(tail -1 "$D/validate/$name.log" | cut -c1-120)"; }
run make_test make test
run dnncheck make dnncheck
run shipcheck make shipcheck
run continuationcheck make continuationcheck
run bench_check python3 bench/bench.py check
DNN_ADMM_CROSS=2 run dnntest_cross python3 tests/dnn_test.py --parts ACDE
DNN_ADMM_WARM_DAMP=0.5:10 DNN_ADMM_WARM_BETA=0.5:10:0.7:1 run dnntest_warm python3 tests/dnn_test.py --parts ACDE
QAP_DNN_WARM_Y=2 DNN_ADMM_CROSS=2 DNN_ADMM_WARM_DAMP=0.5:10 DNN_ADMM_WARM_BETA=0.5:10:0.7:1 \
  run bruteforce_warm python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1
QAP_DNN_WARM_Y=2 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 \
  run bruteforce_wlow python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1
if [ -n "${BASE:-}" ]; then
  echo "== identity_check (BASE=$BASE)"
  BASE=$BASE NEW=$R/latest/dnntest CPU=$CPU bash "$D/identity_check.sh" 2>&1 | tee "$D/validate/identity.log"
fi
