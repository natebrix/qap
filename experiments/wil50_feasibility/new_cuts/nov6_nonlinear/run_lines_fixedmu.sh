#!/bin/bash
# Fixed-multiplier screen of the interior-line box inequalities on node 13
# (README.md section 3): the literature agent's on-top protocol
# (new_cuts/literature/wil50/run_box_ontop.sh), i.e. the grid agent's tri5
# protocol (400 uncut iterations, 5 Polyak rounds of 300 with 1,000 triangles
# a round, target best + 50 uncapped, quiet round of 600) with an external term
# of the 200 most violated interior-line boxes at mu 0.3 (doubled form) added
# to every round's objective, against a matched control run with the same
# binary and settings (bundle D). Both arms on CPU 2, one at a time.
# Predeclared rule: KILL unless the quiet-round certified gain over the control
# is at least +14, the threshold boxes' gain in this screen (literature/
# wil50/runs/n13_ontop_mu0.3_k200.log: 48,728.10 against 48,713.85), a family
# that then failed the pool gate (box_pool/README.md).
set -u
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
Y=/home/user/qap/experiments/wil50_feasibility/new_cuts/grid/runs/n13_tri5.Y.bin
S=21:28,8:41,28:19
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 5 -delta 50 -nocap -freeze 0:600 -x $S"
cd "$D/logs" || exit 1
echo "host $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) amx $(grep -c amx_bf16 /proc/cpuinfo) $(date -u +%FT%TZ)" > host.txt
if ! grep -q '^FINAL' n13_lines_mu0.3_k200.log 2>/dev/null; then
  taskset -c 2 python3 "$D/make_line_ladd.py" $S $Y --mu 0.3 --top 200 --out "$D/logs/n13_lines_L.bin" > n13_lines_mu0.3_k200.mk 2>&1 || { cat n13_lines_mu0.3_k200.mk; exit 1; }
  arg=$(grep -o "use: -Ladd .*" n13_lines_mu0.3_k200.mk | cut -d' ' -f3)
  TIMEFORMAT='PROCESS_CPU user %U sys %S real %R'
  { time taskset -c 2 $C -Ladd "$arg" "$R/qaplib/wil50.dat"; } > n13_lines_mu0.3_k200.log 2>&1
  rm -f "$D/logs/n13_lines_L.bin"
fi
if ! grep -q '^FINAL' n13_tri5_control.log 2>/dev/null; then
  TIMEFORMAT='PROCESS_CPU user %U sys %S real %R'
  { time taskset -c 2 $C "$R/qaplib/wil50.dat"; } > n13_tri5_control.log 2>&1
fi
for f in n13_tri5_control.log n13_lines_mu0.3_k200.log; do
  echo "$f: $(grep '^quiet' $f | grep -o 'lb [0-9.]*') | $(grep '^FINAL' $f) | $(grep PROCESS_CPU $f)"
done
