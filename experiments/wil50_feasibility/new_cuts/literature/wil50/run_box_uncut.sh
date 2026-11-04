#!/bin/bash
# Fixed-multiplier gate of the large-Q 1-box family on the UNCUT DNN (as
# local7/run_local7.sh did for the local-7 cuts): node 13, 1,500 iterations,
# the grid agent's dnntest flags; cuts separated at the uncut iterate
# (grid/runs/n13_uncut.Y.bin, the same flags). Logs runs/n13_unc_*.log.
#   bash run_box_uncut.sh "<mu list>" K
set -u
D=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
MUS=${1:-"1"}; K=${2:-200}
S=21:28,8:41,28:19
T="$R/latest/dnntest -x $S -u 48816 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0"
cd $D/runs
grep -q '^RESULT' n13_unc_ref.log 2>/dev/null || taskset -c 1 $T $R/qaplib/wil50.dat > n13_unc_ref.log 2>&1
[ -s n13_uncut.box.json ] || taskset -c 1 $V -I $D/box_gate.py $R/experiments/wil50_feasibility/new_cuts/grid/runs/n13_uncut --fixes $S --starts 30 --qmax 40 --top $K --out n13_uncut.box.json > n13_uncut.box.log 2>&1
for mu in $MUS; do
  f=n13_unc_box_mu${mu}_k$K.log
  grep -q '^RESULT' $f 2>/dev/null && continue
  arg=$(taskset -c 1 $V -I $D/make_box_ladd.py --json n13_uncut.box.json --mu $mu --top $K --nperm 300 --out n13_L.bin | tee $f.mk | grep -o "use: -Ladd .*" | cut -d' ' -f3)
  [ -n "$arg" ] || { echo "make_box_ladd failed"; cat $f.mk; exit 1; }
  taskset -c 1 $T -Ladd $arg $R/qaplib/wil50.dat > $f 2>&1
  rm -f n13_L.bin
done
for f in n13_unc_ref.log n13_unc_box_*.log; do echo "$f: $(grep -o 'lb=[0-9.]*' $f)"; done
