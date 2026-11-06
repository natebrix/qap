#!/bin/bash
# After the give-up-rule fix of the lazy certificate: tests/dnn_test.py with the
# switch and with the whole bundle on CPU 3, then the all32 run if both pass.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CPU=${CPU:-3}
L=$D/logs/gates
mkdir -p "$L"
cd "$R" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ok=1
for sw in "DNN_CORE_LAZYCERT=1" "DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"; do
  tag=$(echo "$sw" | sed -E 's/DNN_CORE_//g; s/DNN_//g; s/=1//g; s/ /_/g' | tr 'A-Z' 'a-z')
  t0=$(date +%s)
  if env $sw taskset -c $CPU make dnncheck > "$L/dnncheck_${tag}_v2.log" 2>&1; then echo "PASS dnncheck_${tag}_v2 ($(( $(date +%s) - t0 )) s)"; else echo "FAIL dnncheck_${tag}_v2 ($(( $(date +%s) - t0 )) s)"; ok=0; fi
done | tee -a "$L/summary.txt"
grep -q 'FAIL dnncheck_.*_v2' "$L/summary.txt" && { echo "a gate failed: all32 not started" | tee -a "$L/summary.txt"; exit 1; }
exec "$D/run_all32.sh"
