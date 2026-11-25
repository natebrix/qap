#!/bin/bash
# While firings run: every 30 minutes copy their progress and results to snap/ and push it, so a
# reclaimed container loses at most half an hour (copy snap/*.progress.txt back to runs/ to resume).
D=$(cd "$(dirname "$0")"; pwd); cd "$D/../.." || exit 1
while :; do
  cp "$D"/runs/*.progress.txt "$D"/runs/*.sb.txt "$D"/runs/done.txt "$D"/runs/queue.txt "$D/snap/" 2>/dev/null
  git add "$D/snap"
  if ! git diff --cached --quiet && [ "$(git diff --cached | grep -cE '[0-9]{12}')" = 0 ]; then
    git commit -q -m "wil50_sb: snapshot of the firings' progress

Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_014vvMYV3r6FDoHTc2G3VVqC" && git push -q origin claude/cool-babbage-ilnogj
  fi
  pgrep bnbdrv > /dev/null || break
  sleep 1800
done
