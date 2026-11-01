#!/bin/zsh
# freeze_best follow-ups: child 1 with 3 rounds (the best round) and a 2,500 quiet round;
# children 4 (66 short) and 6 (99 short) with 6 rounds and a 2,500 quiet round.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
export VECLIB_MAXIMUM_THREADS=1
C=("$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -u 48816 -T 1 -k 1000 -v -stall 0)
W="$R/qaplib/wil50.dat"
run() { f=$D/$1.log; shift; grep -q '^FINAL' "$f" 2>/dev/null && return; "${C[@]}" "$@" "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
run short3b_d1 -x 39:27,30:48,25:49,7:5  -I 400 -i 300 -r 3 -delta 50 -nocap -freeze 3:2500 &
run short6b_d4 -x 39:27,30:48,25:49,47:5 -I 400 -i 300 -r 6 -delta 50 -nocap -freeze 6:2500 &
run short6b_d6 -x 21:28,8:41,28:19,44:2  -I 400 -i 300 -r 6 -delta 50 -nocap -freeze 6:2500 &
wait
