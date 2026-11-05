#!/bin/bash
# Warm-started children, as in the search (QAP_DNN_WARM_Y=2): two depth-3
# wil50 parents of stronger_bounds/depth4/nodes.txt (clear misses, given up at
# 300 iterations) and six of their depth-4 children, each child started from
# the parent's final Z and conditioned Y (dnntest -K ... -Q b), same flags as
# run_all32.sh. The parent and its children run with the same DNN_RPRSM.
#   usage: run_children.sh BINARY NAME:RPRSM [NAME:RPRSM ...]
# Logs: runs_child/<name>_p<parent>.log.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=$1; shift
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1
mkdir -p "$D/runs_child"
# parent id, parent fixes, children (stronger_bounds/depth4/nodes.txt lines 5-8, 3, 4)
JOBS="13|21:28,8:41,28:19|3:2,44:2,15:2,18:2 10|39:27,30:48,25:49|3:5,47:5"
for cfg in "$@"; do
  name=${cfg%%:*}; rp=${cfg#*:}; [ "$rp" = "$cfg" ] && rp=0
  for job in $JOBS; do
    p=${job%%|*}; rest=${job#*|}; fixes=${rest%%|*}; kids=${rest#*|}
    f=$D/runs_child/${name}_p$p.log
    grep -q '^RESULT node=1 ' "$f" 2>/dev/null && continue
    if [ "$rp" != 0 ]; then export DNN_RPRSM=$rp; else unset DNN_RPRSM; fi
    taskset -c 2 "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 1500 -u 48816 -v \
      -x "$fixes" -K "$kids" -Q b "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
    echo "$name parent $p:"; grep -E '^(CHILD-SPEC|RESULT)' "$f" | sed 's/primal=[^ ]* //'
  done
done
