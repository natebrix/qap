#!/bin/bash
# ./queue.sh ARM...  : the arms one after another (each a 4-worker whole search), after any running one
cd "$(dirname "$0")"
while pgrep -f "bench/parallel.py" > /dev/null; do sleep 5; done
for a in "$@"; do echo "$(date -u +%T) start $a" >> queue.log; ./run.sh $a; echo "$(date -u +%T) end $a: $(cat runs/$a/time.txt)" >> queue.log; done
