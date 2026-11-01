#!/bin/bash
# waits for the clear-miss subtree A/B (search/run_ab2.sh, bench/parallel.py)
# to finish, then runs run_five.sh: at most 4 compute processes at a time
D=$(cd "$(dirname "$0")" && pwd)
while pgrep -f "parallel.py qaplib" > /dev/null || pgrep -x -f "bash run_ab2.sh" > /dev/null || pgrep -x -f "bash ./run_ab2.sh" > /dev/null; do sleep 120; done
echo "A/B done $(date), starting run_five.sh"
bash "$D/run_five.sh"
echo "run_five.sh done $(date)"
