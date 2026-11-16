#!/bin/bash
# core 1 from 04:29: levels 6 and 8 (in case core 0 does not reach them), then the tail again
cd "$(dirname "$0")"
until [ "$(date +%H%M)" -ge 0429 ]; do sleep 10; done
pkill -f "core1d.sh"; pkill -f "queue.sh 1 LA_tail8"; pkill -f "run_dive.sh LA_tail8 229"; pkill -TERM -f "QAP_EST_SEED=229"; pkill -TERM -f "sleep 6000"
P=$(pgrep -f "taskset -c 1 .*A_tail8.state"); [ -n "$P" ] && kill -TERM $P
sleep 2
./queue.sh 1 "LA_off6 226 runs/d2/A_off6.state 6 700 12" "LA_off8 228 runs/d2/A_off8.state 8 700 12" "LA_tail8 249 runs/d2/A_tail8.state 9 6000"
