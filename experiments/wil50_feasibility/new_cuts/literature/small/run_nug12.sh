#!/bin/bash
# nug12 root with the grid on side 0 (as wil50): every relaxation variant, one
# at a time on core 1. Each variant writes its own log with a RESULT line.
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
$V - <<'PY'
import numpy as np, sys
sys.path.insert(0, '.')
from relax_small import read_qaplib, write_qaplib
n, A, B, C = read_qaplib('/home/user/qap/qaplib/nug12.dat')
write_qaplib('inst/nug12_sw.dat', B, A)
PY
for w in ${VARIANTS:-dnn rlt1lp rlt2lp tri lmp1 lmp0 rlt2 lmp01}; do
  [ -f logs/nug12_$w.log ] && grep -q RESULT logs/nug12_$w.log && continue
  taskset -c 1 $V -u relax_small.py inst/nug12_sw.dat --rows 3 --cols 4 --opt 578 --solver ${SOLVER:-SCS} --which $w > logs/nug12_$w.log 2>&1
done
