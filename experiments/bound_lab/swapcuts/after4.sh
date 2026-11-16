#!/bin/sh
# waits for the swapdrv2 run in progress, then runs the job list with swapdrv2: one process at a time
cd "$(dirname "$0")"
while pgrep -x swapdrv2 > /dev/null; do sleep 5; done
exec ./runseq2.sh "$1"
