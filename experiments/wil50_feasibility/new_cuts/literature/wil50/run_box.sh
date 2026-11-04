#!/bin/bash
# Large-Q Jünger-Kaibel 1-box inequalities at a FIXED multiplier on top of the
# regional-triangle phase (the external -Ladd path, no solver changes).
#   ref  : dnncuts triangle protocol (6 rounds of 300 from the 400-iteration
#          give-up point, 1,000 triangles a round, uncapped Polyak delta 50,
#          quiet round 900 from the best round), -dumpY
#   gate : box_gate.py on the ref iterate (13 threshold regions)
#   box  : the same protocol with mu * sum(top-K box cuts) as -Ladd
# One solver process at a time on core 1. Logs runs/<node>_*.log, FINAL = done.
#   bash run_box.sh <node> "<mu list>" [K]
set -u
D=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1
id=$1; MUS=${2:-"1 4"}; K=${3:-50}
spec() { case $1 in d2) echo 39:27,30:48,25:49,26:5;; 13) echo 21:28,8:41,28:19;; 10) echo 39:27,30:48,25:49;; esac; }
S=$(spec $id)
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -v -I 400 -i 300 -r ${ROUNDS:-6} -delta 50 -nocap -freeze 0:${QUIET:-900} -k 1000"
cd $D/runs
if ! grep -q '^FINAL' ${id}_tri.log 2>/dev/null; then
  taskset -c 1 $C -x $S -dumpY ${id}_tri $R/qaplib/wil50.dat > ${id}_tri.log 2>&1
fi
if [ ! -s ${id}_tri.box.json ]; then
  taskset -c 1 $V -I $D/box_gate.py ${id}_tri --fixes $S --starts 30 --qmax 40 --top 200 --out ${id}_tri.box.json > ${id}_tri.box.log 2>&1
fi
for mu in $MUS; do
  f=${id}_tribox_mu${mu}_k${K}.log
  grep -q '^FINAL' $f 2>/dev/null && continue
  arg=$(taskset -c 1 $V -I $D/make_box_ladd.py --json ${id}_tri.box.json --mu $mu --top $K --out ${id}_L_mu${mu}_k${K}.bin | tee $f.mk | grep -o "use: -Ladd .*" | cut -d' ' -f3)
  [ -n "$arg" ] || { echo "make_box_ladd failed"; cat $f.mk; exit 1; }
  taskset -c 1 $C -x $S -Ladd $arg $R/qaplib/wil50.dat > $f 2>&1
  rm -f ${id}_L_mu${mu}_k${K}.bin
done
for f in ${id}_tri.log ${id}_tribox_*.log; do echo "$f: $(grep '^quiet' $f | grep -o 'lb [0-9.]*') | $(grep '^FINAL' $f)"; done
