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
# LAUNCH_TAG=x: a second session's pid file (runs/STAGE.x.pid), e.g. next to a draining one
pidf="$D/runs/$stage${LAUNCH_TAG:+.$LAUNCH_TAG}.pid"
if [ -f "$pidf" ] && kill -0 "$(cat "$pidf")" 2>/dev/null; then
  echo "$stage is already running (pid $(cat "$pidf"))"; exit 1
fi
out="$D/runs/$stage.$(date +%Y%m%d-%H%M%S).out"
cd "$D" || exit 1
setsid nohup python3 "$D/$stage.py" "$@" > "$out" 2>&1 < /dev/null &
echo $! > "$pidf"
echo "started $stage, pid $!, log $out"
