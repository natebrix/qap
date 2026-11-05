#!/bin/bash
# Phase 0 of docs/distributed/plan.md: which CPU to rent (README.md). Run on the
# laptop, signed in to AWS (aws configure sso; aws sso login). It creates, once
# per account and region: a private S3 bucket qap-phase0-<account>-<region>, an
# instance role that can read src/ and write results/ in such buckets and
# nothing else, and a security group with no inbound rules. Then it uploads
# the committed source (git archive of HEAD) and launches one instance per
# family. Each instance builds, benchmarks, uploads to S3 and terminates
# itself; it also shuts itself down after --limit minutes whatever happens.
#
#   infra/phase0/launch.sh [--region eu-north-1] [--families "c7i c7a c7g c8g"]
#                          [--size 4xlarge] [--on-demand] [--serial] [--limit 75]
#                          [--dry-run]
#
#   --serial     one instance at a time (for a vCPU quota below the total)
#   --on-demand  on-demand instead of Spot (about 3-4x the price; same results)
#   --dry-run    print the AWS calls that change something, make none of them
#
# Then: infra/phase0/collect.sh (results and the table), infra/phase0/cleanup.sh.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../.." && pwd)
REGION=eu-north-1
FAMILIES="c7i c7a c7g c8g"
SIZE=4xlarge
MARKET=spot
SERIAL=0
LIMIT=75
DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --region) REGION=$2; shift 2 ;;
    --families) FAMILIES=$2; shift 2 ;;
    --size) SIZE=$2; shift 2 ;;
    --on-demand) MARKET=on-demand; shift ;;
    --serial) SERIAL=1; shift ;;
    --limit) LIMIT=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n 2,21p "$0"; exit 0 ;;
    *) echo "unknown argument $1"; exit 2 ;;
  esac
done
command -v aws > /dev/null || { echo "the AWS CLI is not installed (brew install awscli)"; exit 1; }
export AWS_DEFAULT_REGION=$REGION AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
change() {   # an AWS call that creates or changes something
  if [ "$DRY" = 1 ]; then echo "  [dry-run] aws $*" >&2; echo "dry-run"; else aws "$@"; fi
}
say() { echo "$(date '+%H:%M:%S') $*"; }

ACCT=$(aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
B=qap-phase0-$ACCT-$REGION
RUN=$(date -u +%Y%m%d-%H%M)-$REGION
L=$D/runs/$RUN
[ "$DRY" = 1 ] && L=$(mktemp -d)   # a dry run leaves no run directory
mkdir -p "$L"
say "account $ACCT, region $REGION, run $RUN, families: $FAMILIES ($SIZE, $MARKET)"

# --- quota (vCPUs): Spot "All Standard Spot Instance Requests", on-demand "Running On-Demand Standard"
QCODE=L-34B43A08; [ "$MARKET" = on-demand ] && QCODE=L-1216C47A
Q=$(aws service-quotas get-service-quota --service-code ec2 --quota-code $QCODE --query 'Quota.Value' --output text 2>/dev/null || echo "?")
VCPU=$(aws ec2 describe-instance-types --instance-types "${FAMILIES%% *}.$SIZE" --query 'InstanceTypes[0].VCpuInfo.DefaultVCpus' --output text)
N=$(echo $FAMILIES | wc -w)
NEED=$(( SERIAL == 1 ? VCPU : VCPU * N ))
say "quota $QCODE: $Q vCPUs; this launch needs $NEED ($N x $VCPU$([ "$SERIAL" = 1 ] && echo ', one at a time'))"
if [ "$Q" != "?" ] && [ "${Q%.*}" -lt "$NEED" ]; then
  if [ "$SERIAL" = 0 ] && [ "${Q%.*}" -ge "$VCPU" ]; then
    say "the quota is below $NEED: running one at a time (--serial)"; SERIAL=1
  else
    echo "the $MARKET vCPU quota ($Q) is below one $SIZE instance ($VCPU vCPUs): use a smaller --size, or --on-demand if that quota is higher"; exit 1
  fi
fi

# --- bucket, role, security group (idempotent)
if ! aws s3api head-bucket --bucket "$B" 2> /dev/null; then
  say "creating bucket $B"
  if [ "$REGION" = us-east-1 ]; then change s3api create-bucket --bucket "$B" > /dev/null
  else change s3api create-bucket --bucket "$B" --create-bucket-configuration "LocationConstraint=$REGION" > /dev/null; fi
  change s3api put-public-access-block --bucket "$B" --public-access-block-configuration \
    BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true > /dev/null
fi
ROLE=qap-phase0-instance
if ! aws iam get-role --role-name $ROLE > /dev/null 2>&1; then
  say "creating role $ROLE (EC2 may assume it)"
  change iam create-role --role-name $ROLE --description "qap phase 0 benchmark instances" --assume-role-policy-document \
    '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}' > /dev/null
fi
POLICY="{\"Version\":\"2012-10-17\",\"Statement\":[
 {\"Effect\":\"Allow\",\"Action\":\"s3:GetObject\",\"Resource\":\"arn:aws:s3:::qap-phase0-$ACCT-*/src/*\"},
 {\"Effect\":\"Allow\",\"Action\":\"s3:PutObject\",\"Resource\":\"arn:aws:s3:::qap-phase0-$ACCT-*/results/*\"}]}"
change iam put-role-policy --role-name $ROLE --policy-name qap-phase0-s3 --policy-document "$POLICY" > /dev/null
if ! aws iam get-instance-profile --instance-profile-name $ROLE > /dev/null 2>&1; then
  say "creating instance profile $ROLE"
  change iam create-instance-profile --instance-profile-name $ROLE > /dev/null
  change iam add-role-to-instance-profile --instance-profile-name $ROLE --role-name $ROLE > /dev/null
  [ "$DRY" = 1 ] || sleep 15   # IAM is eventually consistent
fi
VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
[ "$VPC" = None ] && { echo "no default VPC in $REGION (create one: aws ec2 create-default-vpc)"; exit 1; }
SG=$(aws ec2 describe-security-groups --filters Name=group-name,Values=qap-phase0 Name=vpc-id,Values="$VPC" \
  --query 'SecurityGroups[0].GroupId' --output text)
if [ "$SG" = None ]; then
  say "creating security group qap-phase0 in $VPC (no inbound rules)"
  SG=$(change ec2 create-security-group --group-name qap-phase0 --description "qap phase 0: no inbound" \
    --vpc-id "$VPC" --query GroupId --output text)
fi
SUBNETS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values="$VPC" Name=default-for-az,Values=true \
  --query 'Subnets[].SubnetId' --output text)

# --- the source: the committed tree, the parts the build, the gates and the benchmark use
COMMIT=$(git -C "$R" rev-parse --short HEAD)
git -C "$R" diff --quiet HEAD -- latest src infra/phase0 || say "warning: uncommitted changes in latest/, src/ or infra/phase0/ are not uploaded; the instances build $COMMIT"
TAR=$L/qap-$COMMIT.tar.gz
git -C "$R" archive --format=tar.gz -o "$TAR" HEAD Makefile src latest qaplib tests bench coord experiments/dnn_bound \
  experiments/wil50_feasibility/tuning/giveup/all32.txt experiments/wil50_feasibility/speed/r32_noncert/runs infra/phase0
SRC=src/qap-$COMMIT.tar.gz
say "uploading $(du -h "$TAR" | cut -f1) to s3://$B/$SRC"
[ "$DRY" = 1 ] && echo "  [dry-run] aws s3 cp $TAR s3://$B/$SRC" >&2 || aws s3 cp --only-show-errors "$TAR" "s3://$B/$SRC"
rm -f "$TAR"

# --- launch
ami() { aws ssm get-parameter --name "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-$1" \
  --query Parameter.Value --output text; }
AMI_X86=$(ami x86_64); AMI_ARM=$(ami arm64)
echo -e "type\tmarket\tusd_per_hour\tinstance\taz\tlaunched" > "$L/prices.tsv"
launch() {
  local fam=$1 typ=$1.$SIZE ami=$AMI_X86 ud=$L/userdata_$1.sh id="" s
  case $fam in *g|*gd|*gn) ami=$AMI_ARM ;; esac
  cat > "$ud" << EOF
#!/bin/bash
export AWS_DEFAULT_REGION=$REGION QAP_BUCKET=$B QAP_OUT=results/$RUN/$fam
shutdown -h +$LIMIT
mkdir -p /opt/qap && cd /opt/qap
for i in 1 2 3 4 5 6; do aws s3 cp --only-show-errors s3://$B/$SRC src.tgz && break; sleep 10; done
tar xzf src.tgz && bash infra/phase0/bench.sh
shutdown -h now
EOF
  local mo=()
  [ "$MARKET" = spot ] && mo=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
  for s in $SUBNETS; do   # the next zone if one has no capacity
    id=$(change ec2 run-instances --image-id "$ami" --instance-type "$typ" --count 1 --subnet-id "$s" \
      --iam-instance-profile Name=$ROLE --security-group-ids "$SG" --instance-initiated-shutdown-behavior terminate \
      --metadata-options HttpTokens=required,HttpEndpoint=enabled \
      --block-device-mappings 'DeviceName=/dev/xvda,Ebs={VolumeSize=16,VolumeType=gp3,DeleteOnTermination=true}' \
      "${mo[@]}" --user-data "file://$ud" \
      --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=qap-phase0-$fam},{Key=qap-phase0,Value=$RUN}]" \
      --query 'Instances[0].InstanceId' --output text 2> "$L/launch_$fam.err") && break
    say "$typ: not launched in $s ($(tail -1 "$L/launch_$fam.err")); trying the next zone"
    id=""
  done
  [ -n "$id" ] || { say "$typ: could not launch (see $L/launch_$fam.err)"; return 1; }
  local az=? price=?
  if [ "$DRY" != 1 ]; then
    az=$(aws ec2 describe-instances --instance-ids "$id" --query 'Reservations[0].Instances[0].Placement.AvailabilityZone' --output text)
    [ "$MARKET" = spot ] && price=$(aws ec2 describe-spot-price-history --instance-types "$typ" --availability-zone "$az" \
      --product-descriptions Linux/UNIX --start-time "$(date -u +%FT%TZ)" --query 'SpotPriceHistory[0].SpotPrice' --output text)
  fi
  echo -e "$typ\t$MARKET\t$price\t$id\t$az\t$(date -u +%FT%TZ)" >> "$L/prices.tsv"
  say "$typ: $id in $az ($MARKET, \$$price/h); results to s3://$B/results/$RUN/$fam/"
  if [ "$SERIAL" = 1 ] && [ "$DRY" != 1 ]; then
    say "waiting for $id to finish (at most $LIMIT minutes)"
    while [ "$(aws ec2 describe-instances --instance-ids "$id" --query 'Reservations[0].Instances[0].State.Name' --output text)" != terminated ]; do
      sleep 60
    done
  fi
}
for fam in $FAMILIES; do launch "$fam" || true; done
echo "$B" > "$L/bucket.txt"
say "launched. Each instance takes about 30-45 minutes. Then: infra/phase0/collect.sh $RUN"
