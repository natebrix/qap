#!/bin/bash
# commit and push lsla's logs and traces every 20 minutes (SECS, default 40000 s); summ.py refreshes the T2 tables
D=$(cd "$(dirname "$0")"; pwd)
end=$(( $(date +%s) + ${1:-40000} ))
while [ "$(date +%s)" -lt $end ]; do
  ( cd "$D" && for p in tai35b tai40b; do ls -d t2/${p}_* > /dev/null 2>&1 && python3 -I summ.py $p > /dev/null 2>&1; done
    cd "$D/../../.." && git add experiments/branching_lab/lsla &&
    git commit -q -m "branching_lab/lsla: logs in progress $(date -u +%H:%M)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014vvMYSxSbVWLTCUSMMzPcp" && git push -q origin claude/cool-babbage-ilnogj-lsla ) >> "$D/runs/supervise.log" 2>&1
  sleep 1200
done
