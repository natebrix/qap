#!/bin/bash
# second queue (replaces run_nug12.sh after its rlt2lp step): DNN + local
# families by cutting planes, then the LP-level decompositions of RLT2.
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
I=inst/nug12_sw.dat
while [ -n "$WAITPID" ] && kill -0 $WAITPID 2>/dev/null; do sleep 10; done
step() { name=$1; shift; f=logs/nug12_$name.log; grep -q RESULT $f 2>/dev/null && return; taskset -c 1 "$@" > $f 2>&1; }
step rlt2lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which rlt2lp
step tri    $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --solver SCS --eps 1e-5 --which tri
step cp_lmp1 $V -u cp_small.py $I 3 4 lmp1 5 200 1e-5
step cp_lmp0 $V -u cp_small.py $I 3 4 lmp0 5 200 1e-5
step cp_lmp01 $V -u cp_small.py $I 3 4 lmp1+lmp0 5 300 1e-5
step cp_box1 $V -u cp_small.py $I 3 4 box1 5 200 1e-5
step lmp1lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which lmp1lp
step lmp0lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which lmp0lp
step lmp01lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which lmp01lp
step cp_boxall $V -u cp_small.py $I 3 4 boxall 5 300 1e-5
echo QUEUE DONE
