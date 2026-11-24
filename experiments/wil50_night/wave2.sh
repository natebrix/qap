#!/bin/bash
# wil50_night/wave2.sh: the depth-1 line checks (top.sh, base and xl2 at each of the root's 15 children), then
# whatever jobs.txt still holds (queue.sh). Four processes at a time.
D=$(cd "$(dirname "$0")"; pwd); cd "$D"
for k in 1 22 10 11 0 2 3 4 12 13 14 20 21 23 24; do for arm in xl2 base; do
  [ -f "top/${arm}_21-$k.log" ] && continue
  until [ "$(pgrep -c -x bnbdrv)" -lt 4 ]; do sleep 20; done
  nohup ./top.sh $arm 21:$k > /dev/null 2>&1 &
  echo "$(date -u +%T) started top ${arm}_21-$k" >> queue.log
  sleep 5
done; done
wait
echo "$(date -u +%T) top done" >> queue.log
exec ./queue.sh
