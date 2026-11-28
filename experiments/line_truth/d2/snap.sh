#!/bin/bash
# d2/snap.sh: every 30 minutes copy the firings' progress and traces to snap/ and push them, so a lost container
# costs at most that; ends when no firing is left. To resume: copy snap/*.progress.txt to runs/ and run.sh again.
D=$(cd "$(dirname "$0")"; pwd); cd "$D/../../.."
while :; do
  cp "$D"/runs/*.progress.txt "$D"/runs/*.sb.txt "$D"/runs/*.dnn.txt "$D"/runs/done.txt "$D/snap/" 2> /dev/null
  git add "$D/snap" && ! git diff --cached | grep -qE '(^|[^0-9.])[0-9]{12}([^0-9]|$)' && git commit -q -m "line_truth/d2: snapshot of the firings' progress" && git push -q origin claude/cool-babbage-ilnogj
  pgrep -f bnbdrv.d2 > /dev/null || break
  sleep 1800
done
