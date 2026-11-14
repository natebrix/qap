#!/bin/bash
# Phase 1 (README.md, "From GitHub"): let the GitHub workflow .github/workflows/phase1.yml
# launch, collect and stop runs with short-lived AWS credentials (OpenID Connect), so no
# AWS key is stored anywhere. Run once on the laptop, signed in to AWS as for launch.sh and
# to GitHub with the gh CLI (brew install gh; gh auth login). Idempotent: run it again after
# a change to this file.
#
#   infra/phase1/gha_setup.sh [--region eu-north-1] [--dry-run]
#
# AWS:
#   - the OIDC identity provider token.actions.githubusercontent.com (one per account);
#   - a role qap-gha-phase1 that only a job of natebrix/qap running in the GitHub
#     environment aws-phase1 may assume, with the calls launch.sh, collect.sh and
#     cleanup.sh RUN make and nothing else: no IAM writes (launch.sh rewrites the
#     instance role's policy only when it changed, which then needs a launch from the
#     laptop), instances only of the compute families below and only tagged qap-phase1,
#     and only those terminated.
# GitHub (gh):
#   - the environment aws-phase1: you as required reviewer, so every run waits for your
#     approval, and only the cognac branch may deploy to it;
#   - the environment secret AWS_ROLE_ARN (it holds the account number, so it is a secret);
#   - Actions logs and artifacts kept 14 days.
set -euo pipefail
REGION=eu-north-1; DRY=0
REPO=natebrix/qap; ENVNAME=aws-phase1; ROLE=qap-gha-phase1; BRANCH=cognac
FAMILIES='"c6g.*","c6i.*","c7a.*","c7g.*","c7gd.*","c7i.*","c8g.*"'
while [ $# -gt 0 ]; do
  case $1 in
    --region) REGION=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n 2,27p "$0"; exit 0 ;;
    *) echo "unknown argument $1"; exit 2 ;;
  esac
done
command -v aws > /dev/null || { echo "the AWS CLI is not installed (brew install awscli)"; exit 1; }
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
change() { if [ "$DRY" = 1 ]; then echo "  [dry-run] $*" >&2; else "$@"; fi; }
say() { echo "$(date '+%H:%M:%S') $*"; }
ACCT=$(aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }

# --- the identity provider (AWS checks GitHub's certificate itself; the thumbprint is a formality)
PROV=arn:aws:iam::$ACCT:oidc-provider/token.actions.githubusercontent.com
if aws iam get-open-id-connect-provider --open-id-connect-provider-arn "$PROV" > /dev/null 2>&1; then
  say "OIDC provider token.actions.githubusercontent.com: present"
else
  say "creating the OIDC provider token.actions.githubusercontent.com"
  change aws iam create-open-id-connect-provider --url https://token.actions.githubusercontent.com \
    --client-id-list sts.amazonaws.com --thumbprint-list 6938fd4d98bab03faadb97b34396831e3780aea1 > /dev/null
fi

# --- the role: only jobs of $REPO in the environment $ENVNAME
TRUST="{\"Version\":\"2012-10-17\",\"Statement\":[{\"Effect\":\"Allow\",
 \"Principal\":{\"Federated\":\"$PROV\"},\"Action\":\"sts:AssumeRoleWithWebIdentity\",
 \"Condition\":{\"StringEquals\":{\"token.actions.githubusercontent.com:aud\":\"sts.amazonaws.com\",
  \"token.actions.githubusercontent.com:sub\":\"repo:$REPO:environment:$ENVNAME\"}}}]}"
if aws iam get-role --role-name $ROLE > /dev/null 2>&1; then
  say "role $ROLE: present; trust policy rewritten"
  change aws iam update-assume-role-policy --role-name $ROLE --policy-document "$TRUST"
else
  say "creating role $ROLE"
  change aws iam create-role --role-name $ROLE --max-session-duration 3600 \
    --description "qap phase 1 runs from GitHub Actions ($REPO, environment $ENVNAME)" \
    --assume-role-policy-document "$TRUST" > /dev/null
fi

# --- what it may do: the calls of launch.sh, qapctl aws-setup, bk.py, collect.sh, cleanup.sh RUN
B="arn:aws:s3:::qap-coord-$ACCT-*"
INST=qap-phase1-instance
POLICY="{\"Version\":\"2012-10-17\",\"Statement\":[
 {\"Sid\":\"Look\",\"Effect\":\"Allow\",\"Action\":[\"ec2:DescribeInstances\",\"ec2:DescribeInstanceTypes\",\"ec2:DescribeVpcs\",
  \"ec2:DescribeSubnets\",\"ec2:DescribeSecurityGroups\",\"ec2:DescribeSpotPriceHistory\",\"servicequotas:GetServiceQuota\",
  \"dynamodb:ListTables\"],\"Resource\":\"*\"},
 {\"Sid\":\"Image\",\"Effect\":\"Allow\",\"Action\":\"ssm:GetParameter\",
  \"Resource\":\"arn:aws:ssm:*::parameter/aws/service/ami-amazon-linux-latest/*\"},
 {\"Sid\":\"Bucket\",\"Effect\":\"Allow\",\"Action\":\"s3:ListBucket\",\"Resource\":\"$B\"},
 {\"Sid\":\"ReadRuns\",\"Effect\":\"Allow\",\"Action\":\"s3:GetObject\",\"Resource\":[\"$B/runs/*\",\"$B/public/*\"]},
 {\"Sid\":\"Source\",\"Effect\":\"Allow\",\"Action\":\"s3:PutObject\",\"Resource\":\"$B/src/*\"},
 {\"Sid\":\"Queues\",\"Effect\":\"Allow\",\"Action\":[\"sqs:CreateQueue\",\"sqs:GetQueueUrl\",\"sqs:GetQueueAttributes\"],
  \"Resource\":\"arn:aws:sqs:*:$ACCT:qap-units-*\"},
 {\"Sid\":\"Tables\",\"Effect\":\"Allow\",\"Action\":[\"dynamodb:DescribeTable\",\"dynamodb:GetItem\",\"dynamodb:Query\"],
  \"Resource\":[\"arn:aws:dynamodb:*:$ACCT:table/qap-runs\",\"arn:aws:dynamodb:*:$ACCT:table/qap-units\",
   \"arn:aws:dynamodb:*:$ACCT:table/qap-units/index/*\"]},
 {\"Sid\":\"InstanceRoleRead\",\"Effect\":\"Allow\",\"Action\":[\"iam:GetRole\",\"iam:GetRolePolicy\",\"iam:GetInstanceProfile\"],
  \"Resource\":[\"arn:aws:iam::$ACCT:role/$INST\",\"arn:aws:iam::$ACCT:instance-profile/$INST\"]},
 {\"Sid\":\"InstanceRolePass\",\"Effect\":\"Allow\",\"Action\":\"iam:PassRole\",\"Resource\":\"arn:aws:iam::$ACCT:role/$INST\",
  \"Condition\":{\"StringEquals\":{\"iam:PassedToService\":\"ec2.amazonaws.com\"}}},
 {\"Sid\":\"SpotServiceRole\",\"Effect\":\"Allow\",\"Action\":\"iam:CreateServiceLinkedRole\",\"Resource\":\"*\",
  \"Condition\":{\"StringEquals\":{\"iam:AWSServiceName\":\"spot.amazonaws.com\"}}},
 {\"Sid\":\"LaunchTaggedComputeInstances\",\"Effect\":\"Allow\",\"Action\":\"ec2:RunInstances\",
  \"Resource\":\"arn:aws:ec2:*:$ACCT:instance/*\",
  \"Condition\":{\"StringLike\":{\"ec2:InstanceType\":[$FAMILIES]},\"Null\":{\"aws:RequestTag/qap-phase1\":\"false\"}}},
 {\"Sid\":\"LaunchWith\",\"Effect\":\"Allow\",\"Action\":\"ec2:RunInstances\",
  \"Resource\":[\"arn:aws:ec2:*::image/*\",\"arn:aws:ec2:*:$ACCT:subnet/*\",\"arn:aws:ec2:*:$ACCT:network-interface/*\",
   \"arn:aws:ec2:*:$ACCT:security-group/*\",\"arn:aws:ec2:*:$ACCT:volume/*\",\"arn:aws:ec2:*:$ACCT:spot-instances-request/*\"]},
 {\"Sid\":\"TagAtLaunch\",\"Effect\":\"Allow\",\"Action\":\"ec2:CreateTags\",\"Resource\":\"arn:aws:ec2:*:$ACCT:instance/*\",
  \"Condition\":{\"StringEquals\":{\"ec2:CreateAction\":\"RunInstances\"}}},
 {\"Sid\":\"StopOwnInstances\",\"Effect\":\"Allow\",\"Action\":\"ec2:TerminateInstances\",\"Resource\":\"arn:aws:ec2:*:$ACCT:instance/*\",
  \"Condition\":{\"Null\":{\"ec2:ResourceTag/qap-phase1\":\"false\"}}}]}"
python3 -c 'import json, sys; json.loads(sys.argv[1]); json.loads(sys.argv[2])' "$TRUST" "$POLICY" || { echo "bad policy JSON"; exit 1; }
change aws iam put-role-policy --role-name $ROLE --policy-name qap-gha-phase1 --policy-document "$POLICY"
say "role $ROLE: policy qap-gha-phase1 written (no IAM writes; tagged compute instances only)"
ARN=arn:aws:iam::$ACCT:role/$ROLE

# --- GitHub: the environment with you as reviewer, cognac only, the secret, 14-day logs
if ! command -v gh > /dev/null || ! gh auth status > /dev/null 2>&1; then
  echo
  echo "gh is missing or not signed in, so set these by hand on github.com/$REPO:"
  echo "  Settings > Environments > New environment '$ENVNAME': Required reviewers = you;"
  echo "    Deployment branches = Selected branches, add '$BRANCH'"
  echo "  in that environment, Environment secrets > AWS_ROLE_ARN = the ARN of role $ROLE (aws iam get-role --role-name $ROLE)"
  echo "  Settings > Actions > General > Artifact and log retention = 14 days"
  exit 0
fi
ME=$(gh api user --jq .id)
say "GitHub: environment $ENVNAME on $REPO, reviewer $(gh api user --jq .login), branch $BRANCH"
change gh api -X PUT "repos/$REPO/environments/$ENVNAME" --input - > /dev/null << EOF
{"reviewers":[{"type":"User","id":$ME}],"prevent_self_review":false,
 "deployment_branch_policy":{"protected_branches":false,"custom_branch_policies":true}}
EOF
if ! gh api "repos/$REPO/environments/$ENVNAME/deployment-branch-policies" --jq '.branch_policies[].name' 2> /dev/null | grep -qx "$BRANCH"; then
  change gh api -X POST "repos/$REPO/environments/$ENVNAME/deployment-branch-policies" -f name="$BRANCH" -f type=branch > /dev/null
fi
if [ "$DRY" = 1 ]; then echo "  [dry-run] gh secret set AWS_ROLE_ARN --env $ENVNAME --repo $REPO" >&2
else printf '%s' "$ARN" | gh secret set AWS_ROLE_ARN --env "$ENVNAME" --repo "$REPO"; fi
change gh api -X PUT "repos/$REPO/actions/permissions/artifact-and-log-retention" -F days=14 > /dev/null 2>&1 \
  || echo "set the log retention by hand: Settings > Actions > General > Artifact and log retention = 14 days"
say "ready: gh workflow run phase1 --repo $REPO --ref $BRANCH -f action=collect -f run=nug30-live (then approve it)"
