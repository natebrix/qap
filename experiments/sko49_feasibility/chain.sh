#!/bin/bash
# after depth-3 nodes 1-9 finish (10-12 start), stop the depth-3 queue at 12 nodes; when 10-12 finish,
# start the cut runs on node 1 (tri, five) and depth-4 lane B (children 9-14 of node 2)
cd "$(dirname "$0")"
until [ "$(grep -l RESULT depth3/node*.log 2>/dev/null | wc -l)" -ge 9 ]; do sleep 10; done
sleep 30; pkill -f 'xargs -P 3 -L 1'; echo "depth-3 queue stopped at 12 $(date)"
until [ "$(grep -l RESULT depth3/node*.log 2>/dev/null | wc -l)" -ge 12 ]; do sleep 10; done
echo "depth-3 done $(date)"
./cuts/run_cuts.sh tri_1 36:40,4:41,24:10 &
./cuts/run_cuts.sh five_1 36:40,4:41,24:10 -five 500 &
./depth4/run.sh 1 9 14 > depth4/run_b.out 2>&1 &
wait
echo "chain done $(date)"
