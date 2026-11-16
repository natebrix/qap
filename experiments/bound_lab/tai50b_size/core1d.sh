#!/bin/bash
# core 1 after LB_near: the tail below the spine node at depth 8
cd "$(dirname "$0")"
sleep 20
while pgrep -f "queue.sh 1 LB_near" > /dev/null; do sleep 3; done
./queue.sh 1 "LA_tail8 229 runs/d2/A_tail8.state 9 6000"
