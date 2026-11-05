#!/bin/bash
# Phase 0 (README.md): stop and remove what launch.sh made.
#   infra/phase0/cleanup.sh [--region eu-north-1]            terminate every running qap-phase0 instance
#   infra/phase0/cleanup.sh [--region R] --all               also the security group, instance profile and role
#   infra/phase0/cleanup.sh [--region R] --all --bucket      also the bucket and its contents (asks first)
# The instances terminate themselves when they finish or after --limit
# minutes, so the first form is only needed to stop a run early.
set -euo pipefail
REGION=eu-north-1; ALL=0; BUCKET=0
while [ $# -gt 0 ]; do
  case $1 in
    --region) REGION=$2; shift 2 ;;
    --all) ALL=1; shift ;;
    --bucket) BUCKET=1; shift ;;
    *) echo "unknown argument $1"; exit 2 ;;
  esac
done
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
IDS=$(aws ec2 describe-instances --filters Name=tag-key,Values=qap-phase0 Name=instance-state-name,Values=pending,running,stopping,stopped \
  --query 'Reservations[].Instances[].InstanceId' --output text)
if [ -n "$IDS" ]; then echo "terminating $IDS"; aws ec2 terminate-instances --instance-ids $IDS > /dev/null; else echo "no qap-phase0 instances running in $REGION"; fi
[ "$ALL" = 1 ] || exit 0
[ -n "$IDS" ] && { echo "waiting for them to terminate"; aws ec2 wait instance-terminated --instance-ids $IDS || true; }
SG=$(aws ec2 describe-security-groups --filters Name=group-name,Values=qap-phase0 --query 'SecurityGroups[0].GroupId' --output text)
[ "$SG" != None ] && { echo "deleting security group $SG"; aws ec2 delete-security-group --group-id "$SG"; }
ROLE=qap-phase0-instance
if aws iam get-instance-profile --instance-profile-name $ROLE > /dev/null 2>&1; then
  aws iam remove-role-from-instance-profile --instance-profile-name $ROLE --role-name $ROLE 2> /dev/null || true
  aws iam delete-instance-profile --instance-profile-name $ROLE && echo "deleted instance profile $ROLE"
fi
if aws iam get-role --role-name $ROLE > /dev/null 2>&1; then
  aws iam delete-role-policy --role-name $ROLE --policy-name qap-phase0-s3 2> /dev/null || true
  aws iam delete-role --role-name $ROLE && echo "deleted role $ROLE (shared by every region's runs)"
fi
if [ "$BUCKET" = 1 ]; then
  ACCT=$(aws sts get-caller-identity --query Account --output text)
  B=qap-phase0-$ACCT-$REGION
  read -r -p "delete s3://$B and everything in it? [y/N] " a
  [ "$a" = y ] && aws s3 rb "s3://$B" --force
fi
