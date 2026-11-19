#!/bin/bash
# every 15 minutes: rescore one level deeper and commit and push the linetest directory (logs in progress)
D=$(cd "$(dirname "$0")"; pwd)
for i in $(seq 1 16); do
  sleep 900
  ( cd "$D" && python3 -I d9score.py > out/d9score.txt 2>&1; cd "$D/../../.." && git add experiments/branching_lab/linetest &&
    git commit -q -m "branching_lab/linetest: depth-9 logs in progress $(date -u +%H:%M)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014vvMYRPMz8BP8Xc3WigdYm" && git push -q origin claude/cool-babbage-ilnogj-linetest ) > /dev/null 2>&1
done
