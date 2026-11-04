#!/bin/bash
# Sequential job queue on core 3. jobs.txt lines: "<dir> <log> <command ...>"
# (# comments). A job is done when its log has a FINAL line; the log is written
# to <log>.tmp and renamed when the command exits 0. The file is re-read after
# every job, so jobs can be appended while the queue runs. One solver at a time.
D=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=0 DNN_CHOLQR_PASSES=1
while true; do
  next=""
  while read -r dir log cmd; do
    [ -z "$dir" ] && continue; [ "${dir:0:1}" = "#" ] && continue
    grep -q '^FINAL' "$D/$dir/$log" 2>/dev/null && continue
    grep -q '^SKIP' "$D/$dir/$log" 2>/dev/null && continue
    next="$dir|$log|$cmd"; break
  done < "$D/jobs.txt"
  [ -z "$next" ] && { echo "queue empty $(date)"; exit 0; }
  IFS='|' read -r dir log cmd <<< "$next"
  mkdir -p "$D/$dir"; cd "$D/$dir"
  echo "$(date +%H:%M:%S) start $dir/$log: $cmd"
  if taskset -c 3 bash -c "$cmd" > "$log.tmp" 2>&1; then mv "$log.tmp" "$log"; else echo "SKIP failed $(date)" >> "$log.tmp"; mv "$log.tmp" "$log"; fi
  echo "$(date +%H:%M:%S) end $dir/$log: $(grep -m1 '^FINAL' $log)"
  cd "$D"
done
