#!/bin/bash
cd "$(dirname "$0")"; export SOLVER=/home/user/wt-n30/latest/bnbdrv
while pgrep -f "queue12.sh" > /dev/null; do sleep 5; done
echo "$(date -u +%T) start path_g5e" >> queue.log; UB=6126 PAR=g5e ./run.sh path_g5e; echo "$(date -u +%T) end path_g5e: $(cat runs/path_g5e/time.txt)" >> queue.log
