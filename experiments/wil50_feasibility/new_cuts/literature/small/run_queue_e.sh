#!/bin/bash
# queue e: nug12 DNN cutting planes (local RLT2 projections, boxes), then the
# wil50 box test on d2, then the LP-level decompositions.
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
I=inst/nug12_sw.dat
step() { name=$1; shift; f=logs/nug12_$name.log; grep -q RESULT $f 2>/dev/null && return; taskset -c 1 "$@" > $f 2>&1; }
step cp_lmp1 $V -u cp_small.py $I 3 4 lmp1 5 200 1e-5
step cp_lmp0 $V -u cp_small.py $I 3 4 lmp0 5 200 1e-5
step cp_mix  $V -u cp_small.py $I 3 4 mix1+mix0 5 300 1e-5
step cp_box1 $V -u cp_small.py $I 3 4 box1 5 200 1e-5
bash ../wil50/run_box.sh d2 "1 4" 50 > logs/wil50_box_d2.out 2>&1
step cp_lmp01 $V -u cp_small.py $I 3 4 lmp1+lmp0 5 300 1e-5
step lmp1lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which lmp1lp
step cp_boxall $V -u cp_small.py $I 3 4 boxall 5 300 1e-5
echo QUEUE DONE
