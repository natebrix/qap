#!/bin/bash
# usage: run_set.sh CPU NAME FIXES D [-t SECS] key=value ...
# bnbdrv on wil50t from FIXES, bounding nodes of depth <= D only, with the node-10 arm's
# settings (base.par) and the given settings replaced or added. Trace: NAME/trace.txt here.
R=$(cd "$(dirname "$0")/../../../.." && pwd)
Z=$(cd "$(dirname "$0")" && pwd)
CPU=$1; NAME=$2; FIXES=$3; D=$4; shift 4
LIM=""; if [ "$1" = "-t" ]; then LIM="timeout $2"; shift 2; fi
W=$Z/$NAME; mkdir -p "$W"; cp "$Z/base.par" "$W/arm.par"
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}
  if grep -q "^$k = " "$W/arm.par"; then sed -i "s|^$k = .*|$k = $v|" "$W/arm.par"; else echo "$k = $v" >> "$W/arm.par"; fi
done
cd "$W" || exit 1
rm -f trace.txt; date +%H:%M:%S > start.txt
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_DNN_TRACE="$W/trace.txt" $LIM taskset -c "$CPU" "$R/latest/bnbdrv" \
  -p "$W/arm.par" -F "$FIXES" -x 0:4:1:300:1 -s "$R/latest/SYM/wil50t.sym" -u 48816 -d "$D" \
  "$R/qaplib/wil50t.dat" > "$W/solver.log" 2>&1
echo "exit $?" > "$W/done.txt"; date +%H:%M:%S >> "$W/done.txt"
