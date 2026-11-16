#!/bin/bash
# worker.sh ID: take the next line of jobs.txt (under a lock) and run it until the list is empty or
# runs/STOP exists. A job line: SET SEED STATE LEVEL TRIALS. Finished jobs are appended to jobs.done.
H=$(cd "$(dirname "$0")" && pwd); cd $H
J=${JOBS:-jobs.txt}
while [ ! -e runs/STOP ]; do
  job=$( flock 9; j=$(head -n 1 $J); [ -n "$j" ] && sed -i 1d $J; echo "$j" ) 9> runs/jobs.lock
  [ -n "$job" ] || { sleep 30; continue; }
  set -- $job
  ./run_dive.sh $1 $2 $H/runs/strata/$3 $4 $5
  echo "$job | worker $1 $(date -u +%T)" >> jobs.done
done
