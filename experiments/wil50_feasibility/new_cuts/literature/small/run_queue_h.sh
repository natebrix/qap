#!/bin/bash
# after queue g: two more fixed-multiplier box arms on node 13 (uncut)
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
until grep -q "QUEUE G DONE" logs/run_queue_g.out 2>/dev/null; do sleep 20; done
bash ../wil50/run_box_uncut.sh "0.3" 200 > logs/wil50_box_uncut_h1.out 2>&1
bash ../wil50/run_box_uncut.sh "1" 40 > logs/wil50_box_uncut_h2.out 2>&1
echo QUEUE H DONE
