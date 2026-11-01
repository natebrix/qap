#!/bin/zsh
# Follow-up to stepA_d1 (depth-4 child 1 fathomed by the quiet round after 12 rounds):
# shorter protocols on the same child, then the full protocol on child 11 (66 short).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
export VECLIB_MAXIMUM_THREADS=1
C=("$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -u 48816 -T 1 -k 1000 -v -stall 0)
W="$R/qaplib/wil50.dat"
run() { f=$D/$1.log; shift; grep -q '^FINAL' "$f" 2>/dev/null && return; "${C[@]}" "$@" "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
run short4_d1  -x 39:27,30:48,25:49,7:5 -I 400 -i 300 -r 4 -delta 50 -nocap -freeze 4:1500
run short6_d1  -x 39:27,30:48,25:49,7:5 -I 400 -i 300 -r 6 -delta 50 -nocap -freeze 6:1500
run stepA_d11  -x 18:43,11:29,46:44,30:21 -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 12:2000
