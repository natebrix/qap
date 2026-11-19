#!/bin/bash
# replaces the tail of queue5.sh: after the slices, the rest of the opt1 runs (200 s per arm), the had14 cross-check
D=/home/user/qap/experiments/branching_lab/dom_build; cd $D
while pgrep -f "slices.sh tai60b14 tai35b" > /dev/null; do sleep 30; done
python3 run_valid.py logs/valid_opt1.jsonl opt1 200 nug20 had18 nug18 had20 lipa20a lipa20b scr20 tai20b els19 chr18a chr18b chr20a chr20b chr20c esc16c esc16d esc16e esc16f esc16g esc16h esc16i esc16j esc16b > logs/valid_opt1_part3.log 2>&1
BF_MAX=7 python3 xcheck_drops.py logs/xcheck_drops.txt poor had14 tai12b scr12 >> logs/xcheck_drops.log 2>&1
echo QUEUE6 DONE >> logs/valid_opt1_part3.log
