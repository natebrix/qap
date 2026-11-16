#!/bin/bash
# core 0: once the depth-4 spine node (started by spine.sh, which was killed by mistake) is done, the spine to depth 8,
# then the dive sets of the even levels and the tail
cd "$(dirname "$0")"
while pgrep -f "runs/spine_4/in.state" > /dev/null; do sleep 5; done
./spine.sh runs/spine_4/out.state 0 5 8
for d in 5 6 7 8; do python3 splitstate.py runs/spine_$((d-1))/out.state runs/d2/A_off$d.state -offpath x; done
cp runs/spine_8/out.state runs/d2/A_tail8.state
./queue.sh 0 "LA_off4 214 runs/d2/A_off4.state 4 900" "LA_off6 216 runs/d2/A_off6.state 6 900" "LA_off8 218 runs/d2/A_off8.state 8 900" "LA_tail8 219 runs/d2/A_tail8.state 9 3000"
