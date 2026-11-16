#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "queue7.sh" > /dev/null; do sleep 5; done
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
SOLVER=/home/user/wt-n30/latest/bnbdrv UB=6126 PAR=g5 run path_g5
run z060
