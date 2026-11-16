#!/bin/bash
cd "$(dirname "$0")"; export SOLVER=/home/user/wt-n30/latest/bnbdrv
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
run g5e
UB=6126 PAR=g5w144 run path_g5w144
