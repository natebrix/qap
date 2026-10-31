#!/bin/bash
# Regional triangle cuts on wil50's hardest depth-3 node, 34:40,28:22,14:29
# (m = 47), following Astra's protocol and then optimizing the multipliers.
# Four runs at once (one core each), logs in node1/:
#   astra    round 0 = 1500 PRSM iterations uncut, then the 1000 most violated
#            same-region triangles at multiplier 1, 700 iterations, no step
#            (Astra's +39.8 trial)
#   polyak   the same start, then 8 rounds of 700 iterations with Polyak steps
#            toward best + 100, new cuts entering at 0
#   polyak1  as polyak, new cuts entering at multiplier 1
#   control  uncut PRSM for 7100 iterations (1500 + 8 x 700), the equal-budget
#            reference (dnntest -H 0)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd); mkdir -p "$D/node1"
X=34:40,28:22,14:29; W=$R/qaplib/wil50.dat; C=$R/experiments/dnn_cuts/dnncuts
run() { f=$D/node1/$1.log; shift; grep -q '^FINAL\|^RESULT' "$f" 2>/dev/null && return; "$@" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"; }
run astra   $C -grid 5:10 -x $X -u 48816 -T 1 -I 1500 -i 700 -k 1000 -mu0 1 -t 0 -r 1 -v $W &
run polyak  $C -grid 5:10 -x $X -u 48816 -T 1 -I 1500 -i 700 -k 1000 -delta 100 -r 8 -v $W &
run polyak1 $C -grid 5:10 -x $X -u 48816 -T 1 -I 1500 -i 700 -k 1000 -delta 100 -mu0 1 -r 8 -v $W &
run control $R/latest/dnntest -T 1 -d 4 -g 1.0 -G 1.0 -H 0 -i 7100 -u 48816 -v -x $X $W &
wait
