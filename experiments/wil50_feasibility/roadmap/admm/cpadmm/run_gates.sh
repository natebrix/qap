#!/bin/bash
# The merge gates on this tree (README.md, "Side 1"), on cores 1-3 with one thread
# per process: make in latest/ (bnbdrv, dnntest, dnncontinuation), then
#   tests/brute_force.py --jobs 3, tests/resume_test.py, tests/parallel_test.py --jobs 3,
#   tests/coord_test.py, bench/bench.py check, tests/dnn_test.py, tests/ship_test.py,
#   ./latest/dnncontinuation
# Summary (exit code, seconds, last lines) in $OUT.txt, full outputs in $OUT/.
#   bash run_gates.sh [OUT]        (default logs/gates_side1)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
OUT=${1:-$D/logs/gates_side1}
mkdir -p "$OUT"
S=$OUT.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cd "$R" || exit 1
echo "Gates on $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)'), $(date -u +%FT%TZ), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'), taskset -c 1-3, one thread per process, CP switches unset" > "$S"
gate() {   # name command...
  local name=$1; shift
  local t0=$SECONDS
  taskset -c 1-3 "$@" > "$OUT/$name.log" 2>&1
  local rc=$?
  printf '%s: exit %d, %d s, last line: %s\n' "$name" "$rc" $((SECONDS - t0)) "$(grep -v '^ *$' "$OUT/$name.log" | tail -1)" >> "$S"
}
gate make make -C latest bnb dnntest dnncontinuation
gate brute_force python3 tests/brute_force.py --jobs 3
gate resume_test python3 tests/resume_test.py
gate parallel_test python3 tests/parallel_test.py --jobs 3
gate coord_test python3 tests/coord_test.py
gate bench_check python3 bench/bench.py check
gate dnn_test python3 tests/dnn_test.py
gate ship_test python3 tests/ship_test.py
gate dnncontinuation ./latest/dnncontinuation
cat "$S"
