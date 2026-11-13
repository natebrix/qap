#!/bin/bash
# Phase 1 (README.md): serve the coordinator's public status (coord/publish.py,
# public/runs/<run>/ in the run bucket) over HTTPS for brixius-web, without
# making the bucket public. Run on the laptop, signed in to AWS, once per
# account and region; it is idempotent, so running it again only prints the URL.
#
#   infra/phase1/publish_setup.sh [--region eu-north-1] [--dry-run]
#
# It creates:
#   - an origin access control qap-public-<region>: CloudFront signs its requests
#     to the bucket;
#   - a CloudFront distribution whose origin is the bucket with origin path
#     /public, so https://<domain>/runs/<run>/status.json is the object
#     public/runs/<run>/status.json and nothing outside public/ can be reached;
#     GET and HEAD only, HTTPS only, the managed CachingOptimized policy (the
#     objects carry Cache-Control: max-age=60), PriceClass_100;
#   - a bucket policy that lets only that distribution read public/*. The
#     bucket's Block Public Access settings stay as they are: the policy grants
#     a service, not the public.
# A new distribution takes 5 to 15 minutes to deploy. Cost: within CloudFront's
# free tier at a dashboard's traffic (1 TB and 10 million requests a month).
set -euo pipefail
REGION=eu-north-1; DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --region) REGION=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n 2,22p "$0"; exit 0 ;;
    *) echo "unknown argument $1"; exit 2 ;;
  esac
done
command -v aws > /dev/null || { echo "the AWS CLI is not installed (brew install awscli)"; exit 1; }
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
change() { if [ "$DRY" = 1 ]; then echo "  [dry-run] aws $*" >&2; echo "dry-run"; else aws "$@"; fi; }
say() { echo "$(date '+%H:%M:%S') $*"; }

ACCT=$(aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
B=qap-coord-$ACCT-$REGION
aws s3api head-bucket --bucket "$B" 2> /dev/null \
  || { echo "no bucket $B: run infra/phase1/launch.sh (or qapctl aws-setup) first"; exit 1; }
NAME=qap-public-$REGION
COMMENT="qap coordinator public status: s3://$B/public/"

# --- origin access control
OAC=$(aws cloudfront list-origin-access-controls --query "OriginAccessControlList.Items[?Name=='$NAME'].Id | [0]" --output text)
if [ "$OAC" = None ] || [ -z "$OAC" ]; then
  say "creating origin access control $NAME"
  OAC=$(change cloudfront create-origin-access-control --origin-access-control-config \
    "{\"Name\":\"$NAME\",\"Description\":\"qap public status, s3://$B/public/\",\"SigningProtocol\":\"sigv4\",\"SigningBehavior\":\"always\",\"OriginAccessControlOriginType\":\"s3\"}" \
    --query OriginAccessControl.Id --output text)
fi

# --- distribution
read -r DIST DOMAIN <<< "$(aws cloudfront list-distributions \
  --query "DistributionList.Items[?Comment=='$COMMENT'] | [0].[Id,DomainName]" --output text)"
if [ "$DIST" = None ] || [ -z "$DIST" ]; then
  say "creating the CloudFront distribution for s3://$B/public/"
  CFG=$(mktemp)
  cat > "$CFG" << EOF
{
  "CallerReference": "qap-public-$B",
  "Comment": "$COMMENT",
  "Enabled": true,
  "PriceClass": "PriceClass_100",
  "HttpVersion": "http2",
  "Origins": {"Quantity": 1, "Items": [{
    "Id": "qap-public",
    "DomainName": "$B.s3.$REGION.amazonaws.com",
    "OriginPath": "/public",
    "S3OriginConfig": {"OriginAccessIdentity": ""},
    "OriginAccessControlId": "$OAC"
  }]},
  "DefaultCacheBehavior": {
    "TargetOriginId": "qap-public",
    "ViewerProtocolPolicy": "https-only",
    "AllowedMethods": {"Quantity": 2, "Items": ["GET", "HEAD"],
                       "CachedMethods": {"Quantity": 2, "Items": ["GET", "HEAD"]}},
    "CachePolicyId": "658327ea-f89d-4fab-a63d-7e88639e58f6",
    "Compress": true
  }
}
EOF
  read -r DIST DOMAIN <<< "$(change cloudfront create-distribution --distribution-config "file://$CFG" \
    --query 'Distribution.[Id,DomainName]' --output text)"
  [ "$DRY" = 1 ] && { DIST=dry-run; DOMAIN=dry-run.cloudfront.net; }
  rm -f "$CFG"
  # a new account may be refused ("Your account must be verified before you can add new
  # CloudFront resources": an AWS Support case); stop here rather than write a bucket
  # policy for a distribution that does not exist
  [ -n "$DIST" ] && [ "$DIST" != None ] || { echo "the distribution was not created; nothing else changed"; exit 1; }
fi

# --- the bucket policy: that distribution may read public/*, nothing else
POLICY="{\"Version\":\"2012-10-17\",\"Statement\":[{\"Sid\":\"QapPublicStatusViaCloudFront\",\"Effect\":\"Allow\",
 \"Principal\":{\"Service\":\"cloudfront.amazonaws.com\"},\"Action\":\"s3:GetObject\",
 \"Resource\":\"arn:aws:s3:::$B/public/*\",
 \"Condition\":{\"StringEquals\":{\"AWS:SourceArn\":\"arn:aws:cloudfront::$ACCT:distribution/$DIST\"}}}]}"
change s3api put-bucket-policy --bucket "$B" --policy "$POLICY" > /dev/null
STATE=?; [ "$DRY" = 1 ] || STATE=$(aws cloudfront get-distribution --id "$DIST" --query Distribution.Status --output text)
say "distribution $DIST ($STATE; Deployed in 5 to 15 minutes when new)"
echo "public status of run RUN: https://$DOMAIN/runs/RUN/status.json"
echo "for brixius-web, the run's record: \"live\": {\"url\": \"https://$DOMAIN/runs/RUN/\"}"
