#!/bin/bash
# core 1 after LA_off7: dives under the two near-miss B parents alone (31:4, 31:15), the part of the B stratum
# that the estimate expands 45/15; 30 dives
cd "$(dirname "$0")"
while pgrep -f "queue.sh 1 LA_off7" > /dev/null; do sleep 3; done
./queue.sh 1 "LB_near 311 runs/d2/B_near.state 2 4200 30"
