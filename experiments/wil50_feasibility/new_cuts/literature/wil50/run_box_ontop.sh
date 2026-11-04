#!/bin/bash
# Large-Q 1-box cuts at a FIXED multiplier ON TOP of the regional-triangle
# phase, node 13, with exactly the grid agent's protocol of
# grid/runs/n13_tri5.log (-I 400 -i 300 -r 5 -delta 50 -nocap -freeze 0:600,
# 1,000 triangles a round), whose result 48,713.85 is the reference. The box
# cuts are separated at that run's dumped iterate (grid/runs/n13_tri5.Y.bin,
# wil50/n13_tri5.box.json). Logs runs/n13_ontop_mu<MU>_k<K>.log (FINAL = done).
#   bash run_box_ontop.sh "<mu list>" K
set -u
D=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
MUS=${1:-"0.3"}; K=${2:-200}
S=21:28,8:41,28:19
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 5 -delta 50 -nocap -freeze 0:600 -x $S"
cd $D/runs
for mu in $MUS; do
  f=n13_ontop_mu${mu}_k$K.log
  grep -q '^FINAL' $f 2>/dev/null && continue
  arg=$(taskset -c 1 $V -I $D/make_box_ladd.py --json $D/n13_tri5.box.json --mu $mu --top $K --nperm 300 --out n13_Lt.bin | tee $f.mk | grep -o "use: -Ladd .*" | cut -d' ' -f3)
  [ -n "$arg" ] || { echo "make_box_ladd failed"; cat $f.mk; exit 1; }
  taskset -c 1 $C -Ladd $arg $R/qaplib/wil50.dat > $f 2>&1
  rm -f n13_Lt.bin
done
for f in n13_ontop_*.log; do echo "$f: $(grep '^quiet' $f | grep -o 'lb [0-9.]*') | $(grep '^FINAL' $f)"; done
echo "reference (grid/runs/n13_tri5.log): $(grep '^quiet' $R/experiments/wil50_feasibility/new_cuts/grid/runs/n13_tri5.log | grep -o 'lb [0-9.]*')"
