#!/bin/bash
# The default path is unchanged: the same node with a binary built before the
# change ($1) and after it ($2), production-like flags, 150 iterations; the
# trace lines without their CPU column and the RESULT line without cpu must
# be identical. Pinned to core 3, one process at a time.
# usage: identity_check.sh old_dnntest new_dnntest outdir
set -e
R=$(cd "$(dirname "$0")/../../../.." && pwd)
old=$1; new=$2; out=$3
mkdir -p "$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
unset DNN_DUMP_RR DNN_DUMP_RR_AT
for b in old new; do
  bin=$old; [ $b = new ] && bin=$new
  taskset -c 3 "$bin" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 150 -u 48816 -v \
    -x "26:27,8:38,34:28,37:2" "$R/qaplib/wil50.dat" > "$out/ident_$b.log" 2>&1
  grep -E '^(it |RESULT)' "$out/ident_$b.log" | sed -e 's/ cpu [0-9.]*$//' -e 's/ cpu=[0-9.]*//' \
    -e 's/ ms_per_iter=[0-9.]*//' > "$out/ident_$b.trace"
done
if cmp -s "$out/ident_old.trace" "$out/ident_new.trace"; then
  echo "identical: $(wc -l < "$out/ident_new.trace") trace lines"; tail -1 "$out/ident_new.trace"
else
  echo "DIFFERENT"; diff "$out/ident_old.trace" "$out/ident_new.trace" | head
fi
