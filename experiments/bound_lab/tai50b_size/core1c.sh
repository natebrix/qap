#!/bin/bash
# core 1 after LB_near: dives under 31:43 alone (the half-twin; E2 counts it twice, for 31:42)
cd "$(dirname "$0")"
while pgrep -f "queue.sh 1 LB_near" > /dev/null; do sleep 3; done
./queue.sh 1 "LB_43 312 runs/d2/B_43.state 2 3000 20"
