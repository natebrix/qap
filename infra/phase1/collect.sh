#!/bin/bash
# Phase 1 (README.md): where a run stands, and its results once it is done.
#   infra/phase1/collect.sh [RUN]     (default: the latest run under runs/)
# Shows the run's instances, the coordinator's status.json, each host's last
# steps, and results.txt and compare.txt when they exist. Everything it fetches
# goes to runs/RUN/ (each host's out/ under hosts/HOST/). Safe to run at any time.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$(ls -1t "$D/runs" 2> /dev/null | head -1)}
L=$D/runs/$RUN
[ -n "$RUN" ] && [ -f "$L/bucket.txt" ] || { echo "no run ${RUN:-} under $D/runs (launch.sh writes it)"; exit 1; }
B=$(cat "$L/bucket.txt"); REGION=$(cat "$L/region.txt")
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
get compare.txt && { echo; cat "$L/compare.txt"; }
echo
echo "logs: $L/hosts/*/out/ (node.log, coordinator.log, worker.log, gate_*.log)"
