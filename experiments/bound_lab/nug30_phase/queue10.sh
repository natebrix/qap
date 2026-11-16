#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
echo "$(date -u +%T) end path_g5: $(cat runs/path_g5/time.txt)" >> queue.log
echo "$(date -u +%T) start skob" >> queue.log; SOLVER=/home/user/wt-n30/latest/bnbdrv ./run.sh skob; echo "$(date -u +%T) end skob: $(cat runs/skob/time.txt)" >> queue.log
