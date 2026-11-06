#!/bin/bash
# The repository's gates on the committed binary (README section 5), one BLAS thread,
# pinned with taskset. Logs in logs/gates/, one line per gate in logs/gates/summary.txt.
#   SET=single (CPU, default 0): the gates that run one process at a time, so they can
#              run beside a measurement on the other cores:
#     make dnncheck, make shipcheck, make shipbasischeck, make continuationcheck,
#     bench/bench.py check (exact), and with QAP_DNN_WARM_Y=2, DNN_SHIP_BASIS=1 and the
#     production environment: tests/dnn_test.py --parts ACDE (part B compares against
#     the fixed-beta prototype and fails by design under DNN_BETA_SCHED) and
#     tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1; also brute force
#     with the switch alone (the default DNN settings)
#   SET=multi (CPUS, default 0-3): make test (brute force, resume, parallel, coord,
#              regression; its parallel parts run up to 4 processes)
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
L=$D/logs/gates
mkdir -p "$L"
cd "$R" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PROD="DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 DNN_R32_NONCERT=1"
g() { # name cpus command...
  local name=$1 cpus=$2; shift 2
  local t0=$(date +%s)
  if taskset -c "$cpus" "$@" > "$L/$name.log" 2>&1; then echo "PASS $name ($(( $(date +%s) - t0 )) s)"; else echo "FAIL $name ($(( $(date +%s) - t0 )) s): see logs/gates/$name.log"; fi
}
{
echo "# $(date -u +%FT%TZ) gates ($SET) on $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)')"
if [ "${SET:-single}" = single ]; then
  C=${CPU:-0}
  g dnncheck "$C" make dnncheck
  g shipcheck "$C" make shipcheck
  g shipbasischeck "$C" make shipbasischeck
  g continuationcheck "$C" make continuationcheck
  g regression "$C" python3 bench/bench.py check
  g dnn_test_ACDE_prod_ship "$C" env QAP_DNN_WARM_Y=2 DNN_SHIP_BASIS=1 $PROD python3 tests/dnn_test.py --parts ACDE
  g bruteforce_prod_ship "$C" env QAP_DNN_WARM_Y=2 DNN_SHIP_BASIS=1 $PROD python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1
  g bruteforce_ship "$C" env QAP_DNN_WARM_Y=2 DNN_SHIP_BASIS=1 python3 tests/brute_force.py --solver-args "-x 0:20:1:200:1" --jobs 1
else
  g test "${CPUS:-0-3}" make test
fi
echo "# done $(date -u +%T)"
} 2>&1 | tee -a "$L/summary.txt"
