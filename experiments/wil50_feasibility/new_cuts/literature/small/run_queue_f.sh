#!/bin/bash
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
I=inst/nug12_sw.dat
while [ -n "$WAITPID" ] && kill -0 $WAITPID 2>/dev/null; do sleep 10; done
step() { name=$1; shift; f=logs/nug12_$name.log; grep -q RESULT $f 2>/dev/null && return; taskset -c 1 "$@" > $f 2>&1; }
step cp_lmp0 $V -u cp_small.py $I 3 4 lmp0 3 200 1e-5
step cp_mix  $V -u cp_small.py $I 3 4 mix1+mix0 3 300 1e-5
bash ../wil50/run_box_uncut.sh "1 3" 200 > logs/wil50_box_uncut.out 2>&1
step cp_box1 $V -u cp_small.py $I 3 4 box1 3 200 1e-5
step cp_lmp01 $V -u cp_small.py $I 3 4 lmp1+lmp0 3 300 1e-5
echo QUEUE DONE
