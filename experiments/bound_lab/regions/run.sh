#!/bin/bash
# run.sh TAG INSTANCE CUTOFF [dnncuts args...]: one dnncuts run (worktree build), log in logs/TAG.log
# (the verbose per-check lines are filtered to the certified trace: "cp it" and "it" lines every 100)
DC=/home/user/wt-regions/experiments/dnn_cuts/dnncuts
Q=/home/user/qap/qaplib
tag=$1; inst=$2; cut=$3; shift 3
echo "CMD dnncuts -c $cut $* -T 1 $Q/$inst.dat" > logs/$tag.log
$DC -c $cut "$@" -T 1 $Q/$inst.dat 2>&1 | grep -v "^bsched" | awk '
  /^it +[0-9]+ / { if ($2 % 100 == 0) print; next }
  /^cp it +[0-9]+ / { if ($3 % 100 == 0) print; next }
  { print }' >> logs/$tag.log
grep -H "^FINAL" logs/$tag.log
