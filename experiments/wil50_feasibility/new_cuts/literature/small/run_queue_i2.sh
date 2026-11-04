#!/bin/bash
# after queue h: certified fixed-multiplier 1-box test at the nug20 / nug30 roots
# (grid moved to side 0 by swapping the matrices), to see how the family's
# value on top of the DNN scales with n (nug12: +2.3 at mu 2 with 10 cuts,
# full closure with the exact cutting planes).
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
until grep -q "QUEUE H DONE" logs/run_queue_h.out 2>/dev/null; do sleep 20; done
bash ../wil50/run_box_ontop.sh "0.01" 200 > logs/wil50_box_ontop_placebo.out 2>&1
T=/home/user/qap/latest/dnntest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for spec in nug20:4:5:2570 nug30:5:6:6124 nug12:3:4:578; do
  IFS=: read inst r c opt <<< "$spec"
  $V - <<PY
import sys; sys.path.insert(0, '.')
from relax_small import read_qaplib, write_qaplib
n, A, B, C = read_qaplib('/home/user/qap/qaplib/$inst.dat'); write_qaplib('inst/${inst}_sw.dat', B, A)
PY
  [ -s logs/${inst}_sw.Y.bin ] || taskset -c 1 $T -u $opt -i 3000 -H 0 -dumpY logs/${inst}_sw inst/${inst}_sw.dat > logs/${inst}_sw_uncut.log 2>&1
  taskset -c 1 $V -I ../wil50/box_gate.py logs/${inst}_sw --n ${inst#nug} --grid $r:$c --starts 30 --qmax 20 --top 200 --out logs/${inst}_sw.box.json > logs/${inst}_sw.box.log 2>&1
  for K in 20 100; do for mu in 0.5 2; do
    f=logs/${inst}_sw_box_mu${mu}_k$K.log; grep -q RESULT $f 2>/dev/null && continue
    arg=$(taskset -c 1 $V -I ../wil50/make_box_ladd.py --json logs/${inst}_sw.box.json --mu $mu --top $K --nperm 1000 --out logs/L_tmp.bin | tee $f.mk | grep -o "use: -Ladd .*" | cut -d' ' -f3)
    [ -n "$arg" ] && taskset -c 1 $T -u $opt -i 3000 -H 0 -Ladd $arg inst/${inst}_sw.dat > $f 2>&1
  done; done
  rm -f logs/L_tmp.bin
  echo "$inst uncut $(grep -o 'lb=[0-9.]*' logs/${inst}_sw_uncut.log); $(for f in logs/${inst}_sw_box_*.log; do echo -n "$(basename $f .log) $(grep -o 'lb=[0-9.]*' $f); "; done)"
done
echo QUEUE I DONE
