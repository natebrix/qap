#!/bin/bash
# launch.sh stage_a|stage_b [ARGS...]: run the stage detached from this shell
# (its own session, so it outlives the terminal), log in runs/STAGE.<time>.out,
# pid in runs/STAGE.pid. Stop it with: kill -TERM $(cat runs/STAGE.pid)
# (drains: running tasks finish; the same command resumes).
D=$(cd "$(dirname "$0")" && pwd)
stage=$1
shift
case $stage in stage_a|stage_b) ;; *) echo "usage: $0 stage_a|stage_b [ARGS...]"; exit 2;; esac
mkdir -p "$D/runs"
if [ -f "$D/runs/$stage.pid" ] && kill -0 "$(cat "$D/runs/$stage.pid")" 2>/dev/null; then
  echo "$stage is already running (pid $(cat "$D/runs/$stage.pid"))"; exit 1
fi
out="$D/runs/$stage.$(date +%Y%m%d-%H%M%S).out"
cd "$D" || exit 1
setsid nohup python3 "$D/$stage.py" "$@" > "$out" 2>&1 < /dev/null &
echo $! > "$D/runs/$stage.pid"
echo "started $stage, pid $!, log $out"
