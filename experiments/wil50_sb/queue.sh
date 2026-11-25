#!/bin/bash
# Start the next depth-1 firings as cores free up (four at a time), heaviest first by
# their share of wil50-est-261124's estimated rest.
D=$(cd "$(dirname "$0")"; pwd)
for node in "$@"; do
  while [ "$(pgrep -c bnbdrv)" -ge "${JOBS:-4}" ]; do sleep 60; done
  nohup "$D/run.sh" "$node" > /dev/null 2>&1 &
  echo "$node started $(date -u +%FT%TZ)" >> "$D/runs/queue.txt"
  sleep 20
done
wait
