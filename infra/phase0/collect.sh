#!/bin/bash
# Phase 0 (README.md): fetch the results of a run from S3 and write the table.
#   infra/phase0/collect.sh [RUN]     (default: the latest run under runs/)
# Safe to run while instances are still working: it shows their state and
# summarizes whatever has arrived.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
RUN=${1:-$(ls -1 "$D/runs" | tail -1)}
L=$D/runs/$RUN
[ -f "$L/bucket.txt" ] || { echo "no run $RUN under $D/runs"; exit 1; }
B=$(cat "$L/bucket.txt")
REGION=$(echo "$RUN" | cut -d- -f3-)
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
echo "instances of $RUN:"
aws ec2 describe-instances --filters "Name=tag:qap-phase0,Values=$RUN" \
  --query 'Reservations[].Instances[].[Tags[?Key==`Name`]|[0].Value,InstanceId,InstanceType,State.Name]' --output text | sed 's/^/  /'
for fam in $(aws s3 ls "s3://$B/results/$RUN/" | awk '/PRE/ {print $2}' | tr -d /); do
  mkdir -p "$L/$fam"
  aws s3 cp --only-show-errors "s3://$B/results/$RUN/$fam/out.tgz" "$L/$fam/out.tgz" && tar xzf "$L/$fam/out.tgz" -C "$L/$fam" && rm -f "$L/$fam/out.tgz"
  echo "$fam: $(tail -1 "$L/$fam/out/steps.txt")"
done
python3 "$D/summarize.py" table "$L"
echo "To share: git add infra/phase0/runs/$RUN && git commit -m 'phase 0: $RUN' && git push"
