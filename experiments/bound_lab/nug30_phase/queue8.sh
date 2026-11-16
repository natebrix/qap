#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "queue7.sh" > /dev/null; do sleep 5; done
echo "$(date -u +%T) start z060" >> queue.log; ./run.sh z060; echo "$(date -u +%T) end z060: $(cat runs/z060/time.txt)" >> queue.log
