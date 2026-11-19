#!/bin/bash
# the runs in order, one solver process at a time (replaces queue.sh / queue2.sh)
D=/home/user/qap/experiments/branching_lab/dom_build; cd $D
bash identity.sh > logs/identity.log 2>&1
python3 xcheck_drops.py logs/xcheck_drops.txt opt1 tai10b car10pa car10pc car10pf car10pk esc08c chr12a chr12b chr12c had12 had14 nug12 nug14 rou12 scr12 tai12b nug15 scr15 tai15b > logs/xcheck_drops.log 2>&1
bash slices.sh tai35b6 tho40 tai60b17 > logs/slices.log 2>&1
python3 run_valid.py logs/valid_opt1.jsonl opt1 900 $(cat instances_opt1_rest.txt) > logs/valid_opt1_part2.log 2>&1
python3 run_valid.py logs/valid_poor.jsonl poor 600 $(cat instances_poor.txt) > logs/valid_poor.log 2>&1
python3 xcheck_drops.py logs/xcheck_drops.txt poor tai10b car10pa chr12c had14 tai12b scr12 >> logs/xcheck_drops.log 2>&1
bash slices.sh tai60b14 tai35b >> logs/slices.log 2>&1
python3 keyon_mech.py > logs/keyon_mech.log 2>&1
cd /home/user/wt-dombuild
for t in test nodetrace dnncheck continuationcheck shipcheck; do
  ( time make $t ) > $D/logs/repo_checks_$t.log 2>&1; echo "make $t: exit $?" >> $D/logs/repo_checks.log
done
echo QUEUE3 DONE >> $D/logs/repo_checks.log
