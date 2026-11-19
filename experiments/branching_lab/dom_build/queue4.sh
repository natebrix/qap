#!/bin/bash
# after queue3.sh: the rest of the opt1 validity runs, 300 s per arm
D=/home/user/qap/experiments/branching_lab/dom_build; cd $D
until grep -q "QUEUE3 DONE" logs/repo_checks.log 2>/dev/null; do sleep 30; done
python3 run_valid.py logs/valid_opt1.jsonl opt1 300 nug20 chr18a chr18b had18 nug18 els19 chr20a chr20b chr20c had20 lipa20a lipa20b scr20 tai20b esc16c esc16d esc16e esc16f esc16g esc16h esc16i esc16j esc16b rou20 tai20a > logs/valid_opt1_part3.log 2>&1
BF_MAX=7 python3 xcheck_drops.py logs/xcheck_drops.txt poor had14 tai12b scr12 >> logs/xcheck_drops.log 2>&1
echo QUEUE4 DONE >> logs/valid_opt1_part3.log
