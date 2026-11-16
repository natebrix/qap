#!/bin/bash
# core 3 after its queue (LA_off2, LA_off3): level 5 again (its first process finished one dive in 900 s: the second
# dive's first node was still in its phase when the limit hit), then a second process on the tail
cd "$(dirname "$0")"
while pgrep -f "queue.sh 3" > /dev/null; do sleep 5; done
./queue.sh 3 "LA_off5 225 runs/d2/A_off5.state 5 1800" "LA_tail8 239 runs/d2/A_tail8.state 9 6000"
