#!/bin/sh
# waits for the swapdrv run in progress, then runs the job list: one process at a time
cd "$(dirname "$0")"
while pgrep -x swapdrv > /dev/null; do sleep 5; done
exec ./runseq.sh "$1"
