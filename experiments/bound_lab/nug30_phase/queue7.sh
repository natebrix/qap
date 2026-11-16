#!/bin/bash
cd "$(dirname "$0")"; export SOLVER=/home/user/wt-n30/latest/bnbdrv
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
run basea
run g2a
run g4
run g5
run g2s
