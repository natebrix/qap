#!/bin/bash
# Runs the lines of jobs.txt ("INST NAME ARM SPEC [EXTRA]"), in order, keeping at most 4 solver processes
# on the box (jobs.txt may be appended to while this runs). Log: queue.log.
H=$(cd "$(dirname "$0")" && pwd); cd "$H"; k=0
while :; do
  n=$(wc -l < jobs.txt)
  if [ $k -ge $n ]; then [ -f queue.stop ] && break; sleep 20; continue; fi
  if [ "$(pgrep -c -x bnbdrv)" -lt 4 ]; then
    k=$((k+1)); L=$(sed -n "${k}p" jobs.txt); [ -z "$L" ] && continue
    echo "$(date -u +%FT%TZ) start $L" >> queue.log
    ( set -- $L; ./run.sh "$1" "$2" "$3" "$4" "${*:5}"; echo "$(date -u +%FT%TZ) end $L" >> queue.log ) &
    sleep 5
  else sleep 15; fi
done
