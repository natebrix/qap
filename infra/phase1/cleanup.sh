#!/bin/bash
# Phase 1 (README.md): stop and remove what launch.sh made.
#   infra/phase1/cleanup.sh [RUN] [--region eu-north-1]   terminate the run's instances (every qap-phase1 instance without RUN)
#   infra/phase1/cleanup.sh --all                         also the security group, instance profile and role
#   infra/phase1/cleanup.sh --all --storage               also the queues, the tables and the bucket (asks first)
# The instances shut themselves down when the run is finished or after
# --limit minutes, so the first form is only needed to stop one early. A run
# stopped that way can be resumed: launch.sh RUN (after
# `python3 -m coord.qapctl requeue-leased RUN --backend ...` if a worker was
# killed holding units; otherwise their leases run out in 30 minutes).
set -euo pipefail
REGION=eu-north-1; RUN=""; ALL=0; STORAGE=0
while [ $# -gt 0 ]; do
  case $1 in
    --region) REGION=$2; shift 2 ;;
    --all) ALL=1; shift ;;
    --storage) STORAGE=1; shift ;;
    -*) echo "unknown argument $1"; exit 2 ;;
    *) RUN=$1; shift ;;
  esac
done
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
F="Name=tag-key,Values=qap-phase1"; [ -n "$RUN" ] && F="Name=tag:qap-phase1,Values=$RUN"
IDS=$(aws ec2 describe-instances --filters "$F" Name=instance-state-name,Values=pending,running,stopping,stopped \
  --query 'Reservations[].Instances[].InstanceId' --output text)
if [ -n "$IDS" ]; then echo "terminating $IDS"; aws ec2 terminate-instances --instance-ids $IDS > /dev/null
else echo "no qap-phase1 instances${RUN:+ of $RUN} running in $REGION"; fi
[ "$ALL" = 1 ] || exit 0
[ -n "$IDS" ] && { echo "waiting for them to terminate"; aws ec2 wait instance-terminated --instance-ids $IDS || true; }
SG=$(aws ec2 describe-security-groups --filters Name=group-name,Values=qap-phase1 --query 'SecurityGroups[0].GroupId' --output text)
[ "$SG" != None ] && { echo "deleting security group $SG"; aws ec2 delete-security-group --group-id "$SG"; }
ROLE=qap-phase1-instance
if aws iam get-instance-profile --instance-profile-name $ROLE > /dev/null 2>&1; then
  aws iam remove-role-from-instance-profile --instance-profile-name $ROLE --role-name $ROLE 2> /dev/null || true
  aws iam delete-instance-profile --instance-profile-name $ROLE && echo "deleted instance profile $ROLE"
fi
if aws iam get-role --role-name $ROLE > /dev/null 2>&1; then
  aws iam delete-role-policy --role-name $ROLE --policy-name qap-phase1-coord 2> /dev/null || true
  aws iam delete-role --role-name $ROLE && echo "deleted role $ROLE (shared by every region's runs)"
fi
if [ "$STORAGE" = 1 ]; then
  ACCT=$(aws sts get-caller-identity --query Account --output text)
  B=qap-coord-$ACCT-$REGION
  echo "this deletes every coordinator run in $REGION: the ledger (tables qap-runs, qap-units),"
  echo "the queues (qap-units-large, qap-units-small) and s3://$B with all results and logs"
  read -r -p "delete them? [y/N] " a
  [ "$a" = y ] || exit 0
  for q in qap-units-large qap-units-small; do
    U=$(aws sqs get-queue-url --queue-name $q --query QueueUrl --output text 2> /dev/null) && aws sqs delete-queue --queue-url "$U" && echo "deleted queue $q"
  done
  for t in qap-runs qap-units; do
    aws dynamodb delete-table --table-name $t > /dev/null 2>&1 && echo "deleted table $t"
  done
  aws s3 rb "s3://$B" --force && echo "deleted s3://$B"
fi
