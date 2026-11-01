#!/bin/zsh
# Laptop (M2 Pro, Accelerate) runs of research_bound_quality.md section 4 step A
# (inner lag: freeze the multipliers after the last round, 2,000 quiet iterations)
# and of the cut phase's cost on clear misses with the stall rule, from the
# give-up point. Three chains, one core each. Node ids as in run_dual.sh.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
export VECLIB_MAXIMUM_THREADS=1
C=("$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -u 48816 -T 1 -k 1000 -v -stall 0)
W="$R/qaplib/wil50.dat"
run() { f=$D/$1.log; shift; grep -q '^FINAL' "$f" 2>/dev/null && return; "${C[@]}" "$@" "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
# step A at depth 3: step-2 protocol, then freeze 2,000
( run stepA_2  -x 18:43,11:29,46:44 -I 1500 -i 700 -r 8 -delta 100 -nocap -freeze 8:2000
  run clear_10 -x 39:27,30:48,25:49 -I 400 -i 300 -r 12 -delta 50 -nocap -stall 3 -freeze 0:600 ) &
( run stepA_14 -x 47:23,10:44,40:34 -I 1500 -i 700 -r 8 -delta 100 -nocap -freeze 8:2000
  run clear_13 -x 21:28,8:41,28:19 -I 400 -i 300 -r 12 -delta 50 -nocap -stall 3 -freeze 0:600 ) &
# step A at depth 4 (child 1 of depth4/nodes.txt, 62 short), then child 2 (109 short) with the stall rule
( run stepA_d1 -x 39:27,30:48,25:49,7:5 -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 12:2000
  run clear_d2 -x 39:27,30:48,25:49,26:5 -I 400 -i 300 -r 12 -delta 50 -nocap -stall 3 -freeze 0:600 ) &
wait
