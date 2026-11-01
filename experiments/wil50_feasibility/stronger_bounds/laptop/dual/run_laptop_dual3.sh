#!/bin/zsh
# dnncuts with freeze_best (quiet round from the best round's multipliers):
# the step-A protocol again on depth-3 nodes 2 and 14, and the 4-round protocol on depth-4 child 1.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
export VECLIB_MAXIMUM_THREADS=1
C=("$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -u 48816 -T 1 -k 1000 -v -stall 0)
W="$R/qaplib/wil50.dat"
run() { f=$D/$1.log; shift; grep -q '^FINAL' "$f" 2>/dev/null && return; "${C[@]}" "$@" "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
run short4b_d1 -x 39:27,30:48,25:49,7:5 -I 400 -i 300 -r 4 -delta 50 -nocap -freeze 4:1500 &
run stepAb_2   -x 18:43,11:29,46:44 -I 1500 -i 700 -r 8 -delta 100 -nocap -freeze 8:2000 &
run stepAb_14  -x 47:23,10:44,40:34 -I 1500 -i 700 -r 8 -delta 100 -nocap -freeze 8:2000 &
wait
