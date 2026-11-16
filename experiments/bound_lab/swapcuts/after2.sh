#!/bin/sh
# waits for the first queue (after.sh / runseq.sh jobs_c.txt) to end, then runs the second: one process at a time
cd "$(dirname "$0")"
while pgrep -f "jobs_c.txt" > /dev/null; do sleep 10; done
exec ./runseq2.sh "$1"
