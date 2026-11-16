#!/bin/bash
# One worker over stop/jobs.txt (lines: NAME ARM[@ROOTARM] SPEC [EXTRA]; '#' comments), in file order, ONE solver
# process at a time; a job is claimed by mkdir ../runs/.claim2/NAME_ARM; the file is re-read after every job.
H=$(cd "$(dirname "$0")" && pwd); cd "$H"; mkdir -p ../runs/.claim2
while :; do
  got=0
  while read -r NAME ARM SPEC EXTRA; do
    [ -z "${NAME:-}" ] && continue; case $NAME in \#*) continue;; esac
    if mkdir "../runs/.claim2/${NAME}_${ARM%%@*}" 2>/dev/null; then
      echo "$(date -u +%T) start $NAME $ARM" >> seq.log
      ./run2.sh "$NAME" "$ARM" "$SPEC" "$EXTRA" >> seq.log 2>&1
      echo "$(date -u +%T) done  $NAME $ARM" >> seq.log
      got=1; break
    fi
  done < jobs.txt
  if [ $got = 0 ]; then [ -f seq.stop ] && break; sleep 15; fi
done
