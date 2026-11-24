#!/bin/bash
# wil50_night/top15.sh ARM [N]: top.sh ARM at each of the root's 15 children, N at a time (3).
D=$(cd "$(dirname "$0")"; pwd); cd "$D"; arm=$1; n=${2:-3}
for k in 1 22 10 11 0 2 3 4 12 13 14 20 21 23 24; do
  [ -f "top/${arm}_21-$k.log" ] && continue
  until [ "$(pgrep -c -f "top.sh $arm")" -lt "$n" ]; do sleep 20; done
  nohup ./top.sh "$arm" 21:$k > /dev/null 2>&1 &
  sleep 2
done
wait
find top -name "${arm}_*.state" -delete
echo "$(date -u +%T) top15 $arm done" >> queue.log
