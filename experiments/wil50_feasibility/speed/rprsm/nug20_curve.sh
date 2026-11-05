#!/bin/bash
# nug20 root, cold start, production-like flags (degree 3, -L 2, PRSM 1.0/1.0,
# -r 40, one CholQR pass), 400 iterations, no cutoff, no give-up. One run per
# NAME:RPRSM argument (RPRSM empty or 0: switch unset); core 2 only.
#   usage: nug20_curve.sh BINARY NAME:RPRSM [NAME:RPRSM ...]
# Writes logs/nug20_NAME.log; EXTRA adds dnntest flags (e.g. EXTRA=-v for
# the per-check tau lines), ITERS the iteration count.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=$1; shift
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
for cfg in "$@"; do
  NAME=${cfg%%:*}; RP=${cfg#*:}; [ "$RP" = "$cfg" ] && RP=
  if [ -n "$RP" ] && [ "$RP" != 0 ]; then export DNN_RPRSM=$RP; else unset DNN_RPRSM; fi
  taskset -c 2 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 0 -i ${ITERS:-400} -v ${EXTRA:-} \
    "$R/qaplib/nug20.dat" > "$D/logs/nug20_$NAME.log" 2>&1
  echo "$NAME: $(grep RESULT "$D/logs/nug20_$NAME.log")"
  grep '^rprsm bits' "$D/logs/nug20_$NAME.log"
done
