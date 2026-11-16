#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "queue2.sh" > /dev/null; do sleep 5; done
run(){ echo "$(date -u +%T) start $1" >> queue.log; ./run.sh $1; echo "$(date -u +%T) end $1: $(cat runs/$1/time.txt)" >> queue.log; }
run z002
PAR=base run base2
PAR=z003 run z003b
