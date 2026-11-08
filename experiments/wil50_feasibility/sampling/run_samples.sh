#!/bin/bash
# run_samples.sh [F ...]: the forced-root-line samples (README.md, "Root lines"),
# one after the other (default facilities 45, then 25): each bounds the root,
# the 15 depth-1 nodes of facility F in the J1 cells and a seeded uniform sample
# of 200 of the 635 depth-2 nodes, keeps no parent DNN state (--tau 0), and
# drains (exit 4) if its projection exceeds 9 h of wall. The chain stops at the
# first run that does not finish; running it again resumes that run.
# Run detached: bash launch.sh run_samples. To stop: kill -TERM the stage_a.py
# process (it drains), then the chain stops.
D=$(cd "$(dirname "$0")" && pwd)
cd "$D" || exit 1
[ $# -ge 1 ] || set -- 45 25
for F in "$@"; do
  echo "$(date '+%F %T') facility $F: starting"
  python3 "$D/stage_a.py" --run "$D/runs/stageA_f$F" --root-facility "$F" --sample-last 200 --tau 0 --max-wall-h 9
  rc=$?
  echo "$(date '+%F %T') facility $F: stage_a.py exit $rc"
  if [ $rc -ne 0 ]; then
    echo "$(date '+%F %T') stopping the chain (4 = the 9-hour projection guard, 3 = stopped or a failed task)"
    exit $rc
  fi
  python3 "$D/summarize_a.py" "$D/runs/stageA_f$F" --out "$D/results/f$F" > /dev/null
done
echo "$(date '+%F %T') all samples done"
