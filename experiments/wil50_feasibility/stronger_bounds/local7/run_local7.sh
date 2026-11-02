#!/bin/bash
# The local-7 gate on the iterates the triple trial dumped (../triples/<id>.Y.bin,
# the uncut 1,500-iteration iterate), then, if the violations are worth it, the
# fixed-multiplier Lagrangian test through dnntest -Ladd (as in ../triples).
#   bash run_local7.sh <id> <node spec> [subsets per region = 300] [trial threshold = 0.08]
# Logs: <id>.local7.log (the per-region table), <id>.local7.json (summary + the
# 40 most violated cuts), <id>.mu<MU>.log for the trial. One core.
set -u
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
id=$1; spec=$2; K=${3:-300}; thr=${4:-0.08}
Y=$D/../triples/$id
[ -s $Y.Y.bin ] || { echo "no $Y.Y.bin (run ../triples/run_trial.sh $id first, or dnntest -dumpY)"; exit 1; }
cd "$D"
if [ ! -s $id.local7.json ]; then
  python3 local7_separate.py $Y --fixes "$spec" --n 50 --grid 5:10 --subsets $K --seed 1 --dump-top 40 --out $D/$id.local7.json > $id.local7.log 2>&1
fi
tail -3 $id.local7.log
mx=$(python3 -c "import json;print(json.load(open('$D/$id.local7.json'))['overall']['max'])")
if python3 -c "import sys; sys.exit(0 if $mx >= $thr else 1)"; then
  T="$R/latest/dnntest -x $spec -u 48816 -T 1 -d 4 -g 1.0 -G 1.0 -i 1500 -H 0"
  for mu in 1 5 20; do
    f=$id.mu$mu.log; grep -q "^RESULT" $f 2>/dev/null && continue
    arg=$(python3 make_ladd7.py --json $D/$id.local7.json --mu $mu --out $D/$id.L7_$mu.bin | tee $f | grep -o "use: -Ladd .*" | cut -d' ' -f3)
    $T -Ladd $arg $R/qaplib/wil50.dat >> $f 2>&1
  done
  echo "trial: uncut $(grep -o 'lb=[0-9.]*' ../triples/$id.uncut.log | head -1); $(for mu in 1 5 20; do echo -n "mu$mu $(grep -o 'lb=[0-9.]*' $id.mu$mu.log | head -1); "; done)"
else
  echo "max local7 violation $mx < $thr: no trial run"
fi
