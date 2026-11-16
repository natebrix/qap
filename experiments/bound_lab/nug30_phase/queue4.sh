#!/bin/bash
# ./queue4.sh "ARM[:PAR[:UB]]"...  after whatever is running
cd "$(dirname "$0")"
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
for j in "$@"; do IFS=: read a p u <<< "$j"; echo "$(date -u +%T) start $a" >> queue.log; PAR=${p:-$a} UB=${u:-6124} ./run.sh $a; echo "$(date -u +%T) end $a: $(cat runs/$a/time.txt)" >> queue.log; done
