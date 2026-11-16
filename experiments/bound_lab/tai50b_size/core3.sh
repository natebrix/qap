#!/bin/bash
# core 3 after its queue (LA_off2, LA_off3): a second process on the tail
cd "$(dirname "$0")"
while pgrep -f "queue.sh 3" > /dev/null; do sleep 5; done
until [ -s runs/d2/A_tail8.state ]; do sleep 5; done
./queue.sh 3 "LA_tail8 239 runs/d2/A_tail8.state 9 6000"
