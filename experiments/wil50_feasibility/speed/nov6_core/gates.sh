#!/bin/bash
# The repository's gates with this study's switches off (the default path) and on,
# all on CPU 3 (the harnesses' children inherit the affinity). Logs in logs/gates/.
#   make test            brute force, resume, parallel, coord, regression (bench/reference exact)
#   make dnncheck        tests/dnn_test.py, default and with each switch exported
#   brute_force.py       with each switch exported, -x 0:20 (the DNN bound on every node)
#   make shipcheck       state v6 with the parent's DNN state
#   make continuationcheck
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CPU=${CPU:-3}
L=$D/logs/gates
mkdir -p "$L"
cd "$R" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
g() { # name command...
  local name=$1; shift
  local t0=$(date +%s)
  if taskset -c $CPU "$@" > "$L/$name.log" 2>&1; then echo "PASS $name ($(( $(date +%s) - t0 )) s)"; else echo "FAIL $name ($(( $(date +%s) - t0 )) s): see logs/gates/$name.log"; fi
}
{
echo "# $(date -u +%FT%TZ) gates on $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- latest || echo ' (latest/ modified)')"
g test make test
g dnncheck_default make dnncheck
g continuationcheck make continuationcheck
g shipcheck make shipcheck
for sw in "DNN_CORE_LAZYCERT=1" "DNN_CORE_REFRESH_PARTIAL=1" "DNN_CORE_NOSYMV=1 DNN_NOSYMW=1" "DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"; do
  tag=$(echo "$sw" | sed -E 's/DNN_CORE_//g; s/DNN_//g; s/=1//g; s/ /_/g' | tr 'A-Z' 'a-z')
  g "dnncheck_$tag" env $sw make dnncheck
  g "bruteforce_$tag" env $sw python3 tests/brute_force.py --solver-args "-x 0:20" --jobs 1
done
echo "# done $(date -u +%T)"
} 2>&1 | tee -a "$L/summary.txt"
