#!/bin/bash
# Point a running run's worker group at the committed HEAD, from the laptop: the workers the
# group starts from now on build HEAD instead of the commit the run was launched at. The
# coordinating instance and the workers already running keep theirs.
#
#   infra/phase1/roll_workers.sh RUN [--dry-run]
#
# It uploads HEAD's source archive (the same paths launch.sh ships) and writes a new version of
# the group's launch template, a copy of the latest with only the user data's commit changed
# (QAP_COMMIT and the archive's name). The group launches $Latest. Then it lists the group's
# workers with the commit each was started with and prints how to replace the old ones.
#
# For a fix in what a worker runs (Nov 27: coord/worker.py's lease past SQS's 12 hours, for
# wil50-261126). The solver may differ from the run's only where its settings do not: a unit
# resumes a state file only under the same settings hash (latest/state.C), so the script stops
# if latest/settings.C changed since the run's commit.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../.." && pwd)
RUN=${1:-}; DRY=0; [ "${2:-}" = --dry-run ] && DRY=1
[[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]] || { echo "usage: $0 RUN [--dry-run]"; exit 2; }
L=$D/runs/$RUN
REGION=${QAP_REGION:-eu-north-1}; [ -f "$L/region.txt" ] && REGION=$(cat "$L/region.txt")
export AWS_DEFAULT_REGION=$REGION AWS_PAGER=""
ACCT=$(aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
B=qap-coord-$ACCT-$REGION; [ -f "$L/bucket.txt" ] && B=$(cat "$L/bucket.txt")
mask() { sed "s/$ACCT/<account>/g"; }
say() { echo "$(date '+%H:%M:%S') $*" | mask; }
change() { if [ "$DRY" = 1 ]; then echo "  [dry-run] aws $*" | mask >&2; echo "{}"; else aws "$@"; fi; }
ASG=qap-$RUN

COMMIT=$(git -C "$R" rev-parse --short HEAD)
git -C "$R" diff --quiet HEAD -- latest coord bench infra/phase1 \
  || { echo "uncommitted changes in latest/, coord/, bench/ or infra/phase1/: commit them first (the workers build HEAD)"; exit 1; }
UD=$(aws ec2 describe-launch-template-versions --launch-template-name "$ASG" --versions '$Latest' \
  --query 'LaunchTemplateVersions[0].LaunchTemplateData.UserData' --output text) \
  || { echo "no launch template $ASG: is $RUN a run with a worker group?"; exit 1; }
OLD=$(echo "$UD" | base64 -d | sed -n 's/.*QAP_COMMIT=\([0-9a-f]*\).*/\1/p' | head -1)
[ -n "$OLD" ] || { echo "the template's user data names no QAP_COMMIT"; exit 1; }
if [ "$OLD" = "$COMMIT" ]; then say "the template already starts $COMMIT"; else
  git -C "$R" diff --quiet "$OLD" HEAD -- latest/settings.C \
    || { echo "latest/settings.C changed since $OLD: the new workers' settings hash would not match the run's"; exit 1; }
  TAR=$(mktemp -d)/qap-$COMMIT.tar.gz
  git -C "$R" archive --format=tar.gz -o "$TAR" HEAD Makefile latest qaplib QAPSOLN tests bench coord experiments/dnn_bound infra/phase1 \
    'experiments/estimates/*/*.par' experiments/wil50_feasibility/fullrun/fullrun.par
  SRC=src/qap-$COMMIT.tar.gz
  say "uploading $(du -h "$TAR" | cut -f1) to s3://$B/$SRC"
  [ "$DRY" = 1 ] && echo "  [dry-run] aws s3 cp $TAR s3://$B/$SRC" | mask >&2 || aws s3 cp --only-show-errors "$TAR" "s3://$B/$SRC"
  rm -f "$TAR"
  NEW_UD=$(echo "$UD" | base64 -d | python3 -c '
import sys
old, new = sys.argv[1], sys.argv[2]
t = sys.stdin.read()
a, b = f"QAP_COMMIT={old} ", f"src/qap-{old}.tar.gz"
if t.count(a) != 1 or t.count(b) != 1:
    sys.exit(f"the user data does not name {old} once as QAP_COMMIT and once as its archive")
sys.stdout.write(t.replace(a, f"QAP_COMMIT={new} ").replace(b, f"src/qap-{new}.tar.gz"))' "$OLD" "$COMMIT" | base64 | tr -d '\n')
  V=$(change ec2 create-launch-template-version --launch-template-name "$ASG" --source-version '$Latest' \
        --version-description "roll_workers.sh: $OLD to $COMMIT" --launch-template-data "{\"UserData\": \"$NEW_UD\"}" \
      | python3 -c 'import json, sys; print(json.load(sys.stdin).get("LaunchTemplateVersion", {}).get("VersionNumber", "?"))')
  say "launch template $ASG: version $V starts $COMMIT (was $OLD); the group launches it from now on"
fi
echo "workers of $ASG (instance, zone, state, template version):"
aws autoscaling describe-auto-scaling-groups --auto-scaling-group-names "$ASG" \
  --query 'AutoScalingGroups[0].Instances[].[InstanceId,AvailabilityZone,LifecycleState,LaunchTemplate.Version]' --output text | sed 's/^/  /'
cat << EOF
To replace the old workers with new ones (their units come back from their last checkpoints, at
most 5 minutes old, within the 30-minute lease):
  aws autoscaling update-auto-scaling-group --auto-scaling-group-name $ASG --max-size 8 --desired-capacity 8 --region $REGION
then, once the new workers' steps say "worker: 64 slots" (infra/phase1/collect.sh $RUN), each old one:
  aws autoscaling terminate-instance-in-auto-scaling-group --instance-id INSTANCE --should-decrement-desired-capacity --region $REGION
and finally
  aws autoscaling update-auto-scaling-group --auto-scaling-group-name $ASG --max-size 4 --region $REGION
EOF
