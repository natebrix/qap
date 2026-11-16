#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
echo "$(date -u +%T) end z003: $(cat runs/z003/time.txt)" >> queue.log
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
run g1; run sko
UB=6126 PAR=z030 run path_z030
run z060
