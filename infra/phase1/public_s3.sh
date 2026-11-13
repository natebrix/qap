#!/usr/bin/env bash
# Serve the coordinator's public/ prefix straight from S3 while CloudFront is unavailable.
#   public_s3.sh on    public read on public/* (relaxes BlockPublicPolicy and RestrictPublicBuckets;
#                      the two ACL blocks stay on, and object ownership is BucketOwnerEnforced anyway)
#   public_s3.sh off   back to the CloudFront-only policy and all four blocks on
#   public_s3.sh url   print the base URL the run records point at
# Needs AWS_PROFILE set. The CloudFront statement is kept as publish_setup.sh wrote it.
set -euo pipefail
B=${QAP_BUCKET:-qap-coord-033593658551-eu-north-1}
REGION=${QAP_REGION:-eu-north-1}
BASE="https://$B.s3.$REGION.amazonaws.com/public"
cf_stmt() {
  # the existing CloudFront statement, if any, so a later publish_setup.sh run finds what it left
  aws s3api get-bucket-policy --bucket "$B" --query Policy --output text 2>/dev/null \
    | python3 -c 'import json,sys
p=json.load(sys.stdin); s=[x for x in p["Statement"] if x.get("Sid")=="QapPublicStatusViaCloudFront"]
print(json.dumps(s[0]) if s else "")' || true
}
case ${1:-} in
  on)
    CF=$(cf_stmt)
    PUB="{\"Sid\":\"QapPublicStatusDirect\",\"Effect\":\"Allow\",\"Principal\":\"*\",\"Action\":\"s3:GetObject\",\"Resource\":\"arn:aws:s3:::$B/public/*\"}"
    POLICY="{\"Version\":\"2012-10-17\",\"Statement\":[$PUB${CF:+,$CF}]}"
    aws s3api put-public-access-block --bucket "$B" --public-access-block-configuration \
      BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=false,RestrictPublicBuckets=false
    aws s3api put-bucket-policy --bucket "$B" --policy "$POLICY"
    echo "public read on s3://$B/public/*"
    echo "run RUN: $BASE/runs/RUN/status.json ; record url: $BASE/runs/RUN/"
    ;;
  off)
    CF=$(cf_stmt)
    if [ -n "$CF" ]; then
      aws s3api put-bucket-policy --bucket "$B" --policy "{\"Version\":\"2012-10-17\",\"Statement\":[$CF]}"
    else
      aws s3api delete-bucket-policy --bucket "$B"
    fi
    aws s3api put-public-access-block --bucket "$B" --public-access-block-configuration \
      BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true
    echo "public read off; all four public access blocks on"
    ;;
  url) echo "$BASE" ;;
  *) echo "usage: $0 on|off|url" >&2; exit 2 ;;
esac
