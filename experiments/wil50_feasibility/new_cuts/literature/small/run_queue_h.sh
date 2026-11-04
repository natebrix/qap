#!/bin/bash
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
until grep -q "QUEUE G DONE" logs/run_queue_g.out 2>/dev/null; do sleep 20; done
bash ../wil50/run_box_ontop.sh "0.3 1" 200 > logs/wil50_box_ontop.out 2>&1
echo QUEUE H DONE
