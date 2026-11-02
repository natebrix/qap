#!/bin/bash
# Fixed-multiplier triple-consistency trial on one wil50 node. One core.
#   bash run_trial.sh <id> <node spec> [K triples = 60] [n = 50] [ub = 48816]
# Logs: <id>.uncut.log (reference, dumps <id>.Y.bin), <id>.sep.log, <id>.mu<MU>.log,
# <id>.tri_mu<MU>.log (dnncuts step-7 protocol + the term). Resumable by log.
set -u
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
id=$1; spec=$2; K=${3:-60}; n=${4:-50}; ub=${5:-48816}
inst=$R/qaplib/wil50.dat; [ "$n" = 50 ] || inst=$R/qaplib/$6.dat
T="$R/latest/dnntest -x $spec -u $ub -T 1 -d 4 -g 1.0 -G 1.0 -i 1500 -H 0"
cd "$D"
if ! grep -q "^RESULT" $id.uncut.log 2>/dev/null; then
  $T -dumpY $D/$id "$inst" > $id.uncut.log 2>&1
fi
m=$(grep -o "m=[0-9]*" $id.uncut.log | head -1 | cut -d= -f2)
if [ ! -s $id.sep.log ]; then
  python3 pick_triples.py $m $K 1 | while read t; do
    python3 triplet_separate.py $D/$id --fixes "$spec" --facilities "$t" --n $n --seconds 120 | grep -o '"rounded_cut_violation": [0-9.e-]*\|"seconds": [0-9.]*' | tr '\n' ' '; echo " $t"
  done > $id.sep.log 2>&1
fi
for mu in 2 10 50; do
  f=$id.mu$mu.log; grep -q "^RESULT" $f 2>/dev/null && continue
  python3 make_ladd.py --prefix $D/$id --mu $mu --out $D/$id.L$mu.bin > $f 2>&1 || continue
  $T -Ladd $D/$id.L$mu.bin "$inst" >> $f 2>&1
done
# the best fixed mu on top of the regional triangles, step-7 protocol
best=$(for mu in 2 10 50; do echo "$(grep -o 'lb=[0-9.]*' $id.mu$mu.log | head -1 | cut -d= -f2) $mu"; done | sort -g | tail -1 | cut -d' ' -f2)
f=$id.tri_mu$best.log
grep -q "^FINAL" $f 2>/dev/null || $R/experiments/dnn_cuts/dnncuts -x "$spec" -u $ub -grid 5:10 -T 1 -k 1000 -v \
  -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500 -Ladd $D/$id.L$best.bin "$inst" > $f 2>&1
echo "done: $(grep -o 'lb=[0-9.]*' $id.uncut.log | head -1) uncut; $(for mu in 2 10 50; do echo -n "mu$mu $(grep -o 'lb=[0-9.]*' $id.mu$mu.log | head -1); "; done) tri+mu$best $(grep '^FINAL' $f)"
