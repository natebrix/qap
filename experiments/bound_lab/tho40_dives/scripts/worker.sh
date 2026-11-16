#!/bin/bash
# scripts/worker.sh QUEUEFILE: take the first line of QUEUEFILE not yet taken (flock), run it with run.sh, repeat.
# A line "STOP" ends the worker; lines may be appended while workers run.
HERE=$(cd "$(dirname "$0")/.." && pwd); Q=$1
while :; do
  job=$(flock "$Q.lock" python3 - "$Q" <<'PY'
import sys
q=sys.argv[1]; L=open(q).read().split("\n"); out=""
for i,l in enumerate(L):
    if l.strip() and not l.startswith("#"):
        out=l
        if l.strip()!="STOP": L[i]="#taken "+l
        break
open(q,"w").write("\n".join(L)); print(out)
PY
)
  [ -z "$job" ] && { sleep 20; continue; }
  [ "$job" = STOP ] && break
  echo "$(date -u +%T) $job" >> "$HERE/logs/queue.log"
  eval "$HERE/scripts/run.sh $job"
done
