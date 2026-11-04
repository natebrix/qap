#!/bin/bash
# nug12 root with the grid on side 0 (as wil50): relaxation variants, one at a
# time on core 1. Each step writes its own log ending in a RESULT line; a step
# whose log has RESULT is skipped (rerun after a container restart).
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
I=inst/nug12_sw.dat
step() { name=$1; shift; f=logs/nug12_$name.log; grep -q RESULT $f 2>/dev/null && return; taskset -c 1 "$@" > $f 2>&1; }
step rlt1lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which rlt1lp
step rlt2lp $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --which rlt2lp
step tri    $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --solver SCS --eps 1e-5 --which tri
step cp_box1 $V -u cp_small.py $I 3 4 box1 5 200 1e-5
step cp_lmp1 $V -u cp_small.py $I 3 4 lmp1 5 200 1e-5
step cp_lmp0 $V -u cp_small.py $I 3 4 lmp0 5 200 1e-5
step rlt2   $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --solver SCS --eps 1e-4 --which rlt2
step cp_boxall $V -u cp_small.py $I 3 4 boxall 5 300 1e-5
step lmp01  $V -u relax_small.py $I --rows 3 --cols 4 --opt 578 --solver SCS --eps 1e-4 --which lmp01
echo QUEUE DONE
