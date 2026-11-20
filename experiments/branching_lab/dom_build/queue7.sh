#!/bin/bash
# After the container restart of 00:40 UTC Nov 20 (which killed queue6 and the tai35b whole run): the jobs
# that were still to run, started again by the coordinating session. Two processes at once.
D=/home/user/qap/experiments/branching_lab/dom_build; cd $D
( bash slices.sh tai35b >> logs/slices.log 2>&1; echo "TAI35B DONE" >> logs/slices.log ) &
python3 run_valid.py logs/valid_opt1.jsonl opt1 200 nug20 had18 nug18 had20 lipa20a lipa20b scr20 tai20b els19 chr18a chr18b chr20a chr20b chr20c esc16c esc16d esc16e esc16f esc16g esc16h esc16i esc16j esc16b > logs/valid_opt1_part3.log 2>&1
BF_MAX=7 python3 xcheck_drops.py logs/xcheck_drops.txt poor had14 tai12b scr12 >> logs/xcheck_drops.log 2>&1
echo QUEUE6 DONE >> logs/valid_opt1_part3.log
wait
echo QUEUE7 DONE >> logs/slices.log
