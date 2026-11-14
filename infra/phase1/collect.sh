#!/bin/bash
# Phase 1 (README.md): where a run stands, and its results once it is done.
#   infra/phase1/collect.sh [RUN]     (default: the latest run under runs/; a run launched
#                                     elsewhere is fetched into a new runs/RUN/)
# Shows the run's instances, the coordinator's status.json, each host's last
# steps, and results.txt and compare.txt when they exist. Everything it fetches
# goes to runs/RUN/ (each host's out/ under hosts/HOST/). Safe to run at any time.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$(ls -1t "$D/runs" 2> /dev/null | head -1)}
L=$D/runs/$RUN
[ -n "$RUN" ] || { echo "no run under $D/runs: give its name"; exit 1; }
[[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]] || { echo "RUN may have letters, digits, - and _ only"; exit 1; }
# a run launched elsewhere (another clone, or the GitHub workflow) has no directory here
# yet: its region is QAP_REGION's (default eu-north-1), and the store must know the run
REGION=${QAP_REGION:-eu-north-1}
[ -f "$L/region.txt" ] && REGION=$(cat "$L/region.txt")
# bucket.txt stays out of git (the bucket's name holds the account number); without it,
# the bucket is the signed-in account's, as launch.sh names it
if [ -f "$L/bucket.txt" ]; then
  B=$(cat "$L/bucket.txt")
else
  ACCT=$(command aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
  B=qap-coord-$ACCT-$REGION
fi
if [ ! -d "$L" ]; then
  ST=$(python3 "$D/bk.py" "aws:$B?region=$REGION" status "$RUN" 2> /dev/null || echo none)
  [ "$ST" != none ] || { echo "no run $RUN in s3://qap-coord-<account>-$REGION (QAP_REGION for another region)"; exit 1; }
  mkdir -p "$L"
  echo "$REGION" > "$L/region.txt"
  echo "run $RUN ($ST) was launched elsewhere: its records go to $L"
fi
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
echo "instances of $RUN:"
aws ec2 describe-instances --filters "Name=tag:qap-phase1,Values=$RUN" \
  --query 'Reservations[].Instances[].[Tags[?Key==`qap-role`]|[0].Value,InstanceId,InstanceType,InstanceLifecycle||`on-demand`,State.Name,LaunchTime]' \
  --output text | sed 's/^/  /'
get() { aws s3 cp --only-show-errors "s3://$B/runs/$RUN/$1" "$L/$1" 2> /dev/null; }
if get status.json; then
  python3 - "$L/status.json" << 'EOF'
import json, sys, time
s = json.load(open(sys.argv[1]))
age = time.time() - s["updated"]
print(f"status: {s['status']}, {s['wall_s'] / 60:.1f} min, updated {age / 60:.1f} min ago; units {s['units']}, "
      f"queued {s['queue']}; {s['nodes']} nodes ({s['dnn_nodes']} DNN, {s['dnn_fathomed']} DNN-fathomed), "
      f"CPU {s['cpu_s']:.0f} s, incumbent {s['incumbent']}")
EOF
else
  echo "status: no status.json yet (the instance builds and checks the solver first, about 10-15 minutes)"
fi
for h in $(aws s3 ls "s3://$B/runs/$RUN/hosts/" | awk '/PRE/ {print $2}' | tr -d /); do
  mkdir -p "$L/hosts/$h"
  aws s3 cp --only-show-errors "s3://$B/runs/$RUN/hosts/$h/out.tgz" "$L/hosts/$h/out.tgz" \
    && tar xzf "$L/hosts/$h/out.tgz" -C "$L/hosts/$h" && rm -f "$L/hosts/$h/out.tgz"
  echo "host $h, last steps:"
  tail -3 "$L/hosts/$h/out/steps.txt" | cut -c1-300 | sed 's/^/  /'
done
get results.txt && { echo; cat "$L/results.txt"; }
# a finished --trace run: its whole node trace (qapctl nodes; the tree viewer's input) goes
# with its records when it is small enough for git, else it stays in the store
if [ -f "$L/results.txt" ]; then
  NLOG=$(mktemp)
  if (cd "$D/../.." && python3 -m coord.qapctl nodes "$RUN" --backend "aws:$B?region=$REGION" \
        --out "$L/nodes.txt.gz") > "$NLOG" 2>&1; then
    if [ "$(wc -c < "$L/nodes.txt.gz")" -gt 25000000 ]; then
      rm -f "$L/nodes.txt.gz"
      echo "node trace: over 25 MB gzipped, left in the store (qapctl nodes $RUN --out FILE)"
    else
      echo "node trace: $(tail -1 "$NLOG" | sed 's/^[0-9:]* //')"
    fi
  else
    rm -f "$L/nodes.txt.gz"
    grep -q "no frontier/nodes.txt" "$NLOG" || echo "node trace: $(tail -1 "$NLOG")"
  fi
  rm -f "$NLOG"
fi
get compare.txt && { echo; cat "$L/compare.txt"; }
echo
echo "logs: $L/hosts/*/out/ (node.log, coordinator.log, worker.log, gate_*.log)"
