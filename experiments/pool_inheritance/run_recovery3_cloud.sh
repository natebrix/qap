#!/bin/bash
# Fourth attempt at the four continuation (C) runs that lost their executor
# three times (recovery_plan3.json): the same commands, the same source
# (this branch, built from 1f81133), the committed parent pool, on the Fable
# cloud box (4 x Xeon @ 2.1 GHz, AMX-capable but unused here; scipy-openblas32
# 0.3.34 against 0.3.30 in the original). One process per core. Resumable at
# the run level: a run whose log has a FINAL line is skipped; an interrupted
# run restarts from its start (each is at most 5,500 iterations).
# Reproduction check first (s11_C_repro.log): same pool load, same fathom
# iteration 160, certificate 48,819.83 against 48,819.77.
D=$(cd "$(dirname "$0")" && pwd); W=$(cd "$D/../.." && pwd)
OUT=$D/results/recovery3_cloud; mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
POOL=$D/results/runs/parent.pool
one() {
  name=$1; fix=$2; f=$OUT/${name}_C.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  cd "$W" || exit 1
  # whole-process user + system CPU from the shell's timer (no /usr/bin/time here)
  { TIMEFORMAT='PROCESS_CPU user %U sys %S real %R'; time experiments/dnn_cuts/dnncuts -grid 5:10 -side 0 -u 48816 -T 1 -v -delta 50 -nocap \
    -x 39:27,30:48,25:49,$fix -I 400 -r 12 -i 300 -k 1000 -freeze 0:1500 -loadpool "$POOL" qaplib/wil50.dat; } \
    > "$f.tmp" 2>&1 && grep -q '^FINAL' "$f.tmp" && mv "$f.tmp" "$f"
}
export -f one; export OUT POOL W
printf '%s\n' "d4 47:5" "s5 5:5" "d2 26:5" "d1 7:5" | xargs -P 4 -L 1 bash -c 'one "$@"' _
