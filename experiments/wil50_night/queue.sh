#!/bin/bash
# wil50_night/queue.sh: keeps four solver processes running from jobs.txt (lines: NAME SPEC ARM), in order.
D=$(cd "$(dirname "$0")"; pwd); cd "$D"
while read -r name spec arm; do
  [ -n "$name" ] || continue
  [ -f "runs/${name}_$arm/log.txt" ] && continue
  until [ "$(pgrep -c -x bnbdrv)" -lt 4 ]; do sleep 20; done
  nohup ./run.sh "$name" "$spec" "$arm" > /dev/null 2>&1 &
  echo "$(date -u +%T) started ${name}_$arm" >> queue.log
  sleep 5
done < jobs.txt
wait
echo "$(date -u +%T) queue done" >> queue.log
