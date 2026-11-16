#!/bin/bash
# after whatever is running; W: the worktree build with cut.giveup_depth (same source otherwise)
cd "$(dirname "$0")"; W=/home/user/wt-n30/latest/bnbdrv
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
SOLVER=$W run g2
PAR=base run base2
run sko
[ -f queue6.sh ] && exec ./queue6.sh
