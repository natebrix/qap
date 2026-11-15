#!/bin/bash
# Two workers over jobs.txt (lines: NAME ARM SPEC; '#' comments), in file order; a job is
# claimed by mkdir runs/.claim/NAME_ARM. Re-reads jobs.txt after every job, so lines can be
# appended while it runs. Stops when a pass finds nothing to claim and the file "queue.stop" exists,
# else sleeps and retries.
H=$(cd "$(dirname "$0")" && pwd); cd "$H"; mkdir -p runs/.claim
worker() {
  while :; do
    got=0
    while read -r NAME ARM SPEC EXTRA; do
      [ -z "${NAME:-}" ] && continue; case $NAME in \#*) continue;; esac
      if mkdir "runs/.claim/${NAME}_$ARM" 2>/dev/null; then
        echo "$(date -u +%T) w$1 start $NAME $ARM" >> queue.log
        ./run.sh "$NAME" "$ARM" "$SPEC" "$EXTRA" >> queue.log 2>&1
        echo "$(date -u +%T) w$1 done  $NAME $ARM" >> queue.log
        got=1; break
      fi
    done < jobs.txt
    if [ $got = 0 ]; then [ -f queue.stop ] && break; sleep 20; fi
  done
}
worker 1 & worker 2 & wait
