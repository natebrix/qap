#!/bin/bash
# After a container restart: restart the job queue if it is not running.
D=$(cd "$(dirname "$0")" && pwd)
if ps -eo cmd | grep -v grep | grep -q "bash ./runq.sh"; then echo "queue running"; else
  (cd "$D" && setsid nohup bash ./runq.sh >> runq.log 2>&1 &); echo "queue restarted"; fi
