#!/bin/bash
# Phase 1 of docs/distributed/plan.md (README.md): the coordinator and a worker
# on one EC2 instance, on S3, SQS and DynamoDB. Run on the laptop, signed in to
# AWS (aws sso login), from a clone with boto3 for python3 (pip3 install boto3).
#
#   infra/phase1/launch.sh RUN [--recipe nug30] [--no-ub] [--type c7g.4xlarge]
#                          [--region eu-north-1] [--on-demand] [--worker-only]
#                          [--no-reference] [--gates "bruteforce ..."|none]
#                          [--limit 120] [--node-env "K=V ..."] [--dry-run]
#
#   RUN             a name for the run: letters, digits, - and _ (e.g. nug30-a)
#   --recipe R      what to solve: infra/phase1/recipes/R.sh (committed). Needed
#                   only for a new run; launching an existing RUN resumes it
#   --no-ub         start without the recipe's incumbent: the run must find the optimum
#   --type T        the instance type (default c7g.4xlarge: 16 Graviton3 cores, 16 vCPUs)
#   --on-demand     on-demand instead of Spot (about 3-4x the price; no interruptions)
#   --worker-only   a worker for RUN, which another instance coordinates
#   --no-reference  skip the bench/parallel.py run and the comparison at the end
#   --private       publish no public status for a new run (coord/publish.py;
#                   publish_setup.sh serves it to brixius-web)
#   --gates G       the make targets checked before any unit runs (default in node.sh)
#   --limit M       the instance shuts itself down after M minutes whatever happens
#   --node-env E    variables exported on the instance for the build, the gates and the
#                   solvers (e.g. "OPENBLAS_CORETYPE=ARMV8"); the recipe's ENV is for units only
#   --dry-run       print the AWS calls that change something, make none of them
#
# It creates, once per account and region, and keeps for later runs: the
# coordinator's storage (`qapctl aws-setup`: a private bucket
# qap-coord-<account>-<region>, the SQS queues qap-units-large and -small, the
# DynamoDB tables qap-runs and qap-units), an instance role that may use those
# and nothing else, and a security group with no inbound rules. Then it uploads
# the committed source and launches the instance (node.sh does the rest).
# Then: infra/phase1/collect.sh RUN (progress and results), cleanup.sh.
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../.." && pwd)
RUN=""; RECIPE=""; NO_UB=0; TYPE=c7g.4xlarge; REGION=eu-north-1; MARKET=spot; ROLE=all
REFERENCE=1; GATES=""; LIMIT=120; DRY=0; VIS=public; NODE_ENV=""
while [ $# -gt 0 ]; do
  case $1 in
    --recipe) RECIPE=$2; shift 2 ;;
    --no-ub) NO_UB=1; shift ;;
    --type) TYPE=$2; shift 2 ;;
    --region) REGION=$2; shift 2 ;;
    --on-demand) MARKET=on-demand; shift ;;
    --worker-only) ROLE=worker; shift ;;
    --no-reference) REFERENCE=0; shift ;;
    --private) VIS=private; shift ;;
    --gates) GATES=$2; shift 2 ;;
    --limit) LIMIT=$2; shift 2 ;;
    --node-env) NODE_ENV=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n 2,31p "$0"; exit 0 ;;
    -*) echo "unknown argument $1"; exit 2 ;;
    *) [ -z "$RUN" ] || { echo "one RUN only ($RUN, $1)"; exit 2; }; RUN=$1; shift ;;
  esac
done
[ -n "$RUN" ] || { sed -n 2,31p "$0"; exit 2; }
[[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]] || { echo "RUN may have letters, digits, - and _ only"; exit 2; }
[[ $LIMIT =~ ^[0-9]+$ ]] || { echo "--limit is in whole minutes"; exit 2; }
command -v aws > /dev/null || { echo "the AWS CLI is not installed (brew install awscli)"; exit 1; }
python3 -c "import boto3" 2> /dev/null || { echo "python3 has no boto3: pip3 install boto3"; exit 1; }
export AWS_DEFAULT_REGION=$REGION AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
change() {   # an AWS call that creates or changes something
  if [ "$DRY" = 1 ]; then echo "  [dry-run] aws $*" >&2; echo "dry-run"; else aws "$@"; fi
}
say() { echo "$(date '+%H:%M:%S') $*"; }

ACCT=$(aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
B=qap-coord-$ACCT-$REGION
BK="aws:$B?region=$REGION${QAP_AWS_ENDPOINT:+&endpoint=$QAP_AWS_ENDPOINT}"   # QAP_AWS_ENDPOINT: tests only (moto)
bk() { python3 "$D/bk.py" "$BK" "$@"; }
say "account $ACCT, region $REGION, run $RUN ($ROLE), $TYPE $MARKET"

# --- quota (vCPUs): Spot "All Standard Spot Instance Requests", on-demand "Running On-Demand Standard"
QCODE=L-34B43A08; [ "$MARKET" = on-demand ] && QCODE=L-1216C47A
Q=$(aws service-quotas get-service-quota --service-code ec2 --quota-code $QCODE --query 'Quota.Value' --output text 2>/dev/null || echo "?")
VCPU=$(aws ec2 describe-instance-types --instance-types "$TYPE" --query 'InstanceTypes[0].VCpuInfo.DefaultVCpus' --output text)
say "quota $QCODE: $Q vCPUs; $TYPE needs $VCPU"
if [ "$Q" != "?" ] && [ "${Q%.*}" -lt "$VCPU" ]; then
  echo "the $MARKET vCPU quota in $REGION ($Q) is below one $TYPE ($VCPU vCPUs). Check the request:"
  echo "  aws service-quotas list-requested-service-quota-change-history --service-code ec2 --region $REGION"
  echo "or try another --region, --type, or --on-demand if that quota is higher"
  [ "$DRY" = 1 ] || exit 1
  say "dry run: going on to show the rest (a real launch stops here)"
fi

# --- the coordinator's storage (idempotent; the same definitions as the moto tests)
if [ "$DRY" = 1 ]; then echo "  [dry-run] python3 -m coord.qapctl aws-setup --backend '$BK'" >&2
else (cd "$R" && python3 -m coord.qapctl aws-setup --backend "$BK" > /dev/null); say "storage ready: s3://$B, queues qap-units-*, tables qap-runs and qap-units"; fi

# --- the run: new, resumed, or finished; at most one coordinating instance
ST=$(bk status "$RUN" 2> /dev/null || echo none)
live() {   # instances of this run in a role that are pending or running
  aws ec2 describe-instances --filters "Name=tag:qap-phase1,Values=$RUN" "Name=tag:qap-role,Values=$1" \
    Name=instance-state-name,Values=pending,running --query 'Reservations[].Instances[].InstanceId' --output text
}
COORD_LIVE=$(live all)
OTHERS=$(aws ec2 describe-instances --filters Name=tag:qap-role,Values=all Name=instance-state-name,Values=pending,running \
  --query "Reservations[].Instances[].Tags[?Key=='qap-phase1'].Value[]" --output text | tr -s '[:space:]' '\n' \
  | grep -vx -e "$RUN" -e '' -e None | sort -u | tr '\n' ' ' || true)
[ -z "${OTHERS// /}" ] || say "note: run(s) $OTHERS also running in $REGION; runs share the queues, so each is slower while they overlap"
if [ "$ROLE" = all ]; then
  [ -z "$COORD_LIVE" ] || { echo "run $RUN already has a coordinating instance ($COORD_LIVE); add capacity with --worker-only"; exit 1; }
  case $ST in
    finished) echo "run $RUN is finished: infra/phase1/collect.sh $RUN"; exit 1 ;;
    none) [ -n "$RECIPE" ] || { echo "run $RUN is new: give it a --recipe ($(ls "$D/recipes" | sed 's/\.sh$//' | tr '\n' ' '))"; exit 1; }
          RF=infra/phase1/recipes/$RECIPE.sh
          [ -f "$R/$RF" ] || { echo "no recipe $RF"; exit 1; }
          git -C "$R" ls-files --error-unmatch "$RF" > /dev/null 2>&1 && git -C "$R" diff --quiet HEAD -- "$RF" \
            || { echo "$RF is not committed: the instance builds from HEAD, so commit it first"; exit 1; }
          say "new run $RUN: recipes/$RECIPE.sh$([ "$NO_UB" = 1 ] && echo ', no incumbent')" ;;
    *) say "run $RUN exists ($ST): resuming it$([ -n "$RECIPE" ] && echo '; --recipe is ignored')" ;;
  esac
else
  [ "$ST" != none ] || [ -n "$COORD_LIVE" ] || { echo "no run $RUN and no instance coordinating it: launch one without --worker-only"; exit 1; }
  [ "$ST" != finished ] || { echo "run $RUN is finished"; exit 1; }
fi

# --- instance role, profile and security group (idempotent)
ROLE_NAME=qap-phase1-instance
if ! aws iam get-role --role-name $ROLE_NAME > /dev/null 2>&1; then
  say "creating role $ROLE_NAME (EC2 may assume it)"
  change iam create-role --role-name $ROLE_NAME --description "qap phase 1 coordinator and worker instances" --assume-role-policy-document \
    '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}' > /dev/null
fi
# the calls coord/backend.py makes once the storage exists, on the qap-* names only
POLICY="{\"Version\":\"2012-10-17\",\"Statement\":[
 {\"Effect\":\"Allow\",\"Action\":\"s3:ListBucket\",\"Resource\":\"arn:aws:s3:::qap-coord-$ACCT-*\"},
 {\"Effect\":\"Allow\",\"Action\":\"s3:GetObject\",\"Resource\":\"arn:aws:s3:::qap-coord-$ACCT-*/src/*\"},
 {\"Effect\":\"Allow\",\"Action\":[\"s3:GetObject\",\"s3:PutObject\",\"s3:DeleteObject\",\"s3:AbortMultipartUpload\"],
  \"Resource\":\"arn:aws:s3:::qap-coord-$ACCT-*/runs/*\"},
 {\"Effect\":\"Allow\",\"Action\":\"s3:PutObject\",\"Resource\":\"arn:aws:s3:::qap-coord-$ACCT-*/public/runs/*\"},
 {\"Effect\":\"Allow\",\"Action\":[\"sqs:ReceiveMessage\",\"sqs:SendMessage\",\"sqs:DeleteMessage\",\"sqs:ChangeMessageVisibility\",
  \"sqs:GetQueueUrl\",\"sqs:GetQueueAttributes\"],\"Resource\":\"arn:aws:sqs:*:$ACCT:qap-units-*\"},
 {\"Effect\":\"Allow\",\"Action\":[\"dynamodb:GetItem\",\"dynamodb:PutItem\",\"dynamodb:UpdateItem\",\"dynamodb:Query\"],
  \"Resource\":[\"arn:aws:dynamodb:*:$ACCT:table/qap-runs\",\"arn:aws:dynamodb:*:$ACCT:table/qap-units\",
               \"arn:aws:dynamodb:*:$ACCT:table/qap-units/index/*\"]}]}"
change iam put-role-policy --role-name $ROLE_NAME --policy-name qap-phase1-coord --policy-document "$POLICY" > /dev/null
if ! aws iam get-instance-profile --instance-profile-name $ROLE_NAME > /dev/null 2>&1; then
  say "creating instance profile $ROLE_NAME"
  change iam create-instance-profile --instance-profile-name $ROLE_NAME > /dev/null
  change iam add-role-to-instance-profile --instance-profile-name $ROLE_NAME --role-name $ROLE_NAME > /dev/null
  [ "$DRY" = 1 ] || sleep 15   # IAM is eventually consistent
fi
VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
[ "$VPC" = None ] && { echo "no default VPC in $REGION (create one: aws ec2 create-default-vpc)"; exit 1; }
SG=$(aws ec2 describe-security-groups --filters Name=group-name,Values=qap-phase1 Name=vpc-id,Values="$VPC" \
  --query 'SecurityGroups[0].GroupId' --output text)
if [ "$SG" = None ]; then
  say "creating security group qap-phase1 in $VPC (no inbound rules)"
  SG=$(change ec2 create-security-group --group-name qap-phase1 --description "qap phase 1: no inbound" \
    --vpc-id "$VPC" --query GroupId --output text)
fi
SUBNETS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values="$VPC" Name=default-for-az,Values=true \
  --query 'Subnets[].SubnetId' --output text)

# --- the source: the committed tree, the parts the build, the gates and the run use
COMMIT=$(git -C "$R" rev-parse --short HEAD)
git -C "$R" diff --quiet HEAD -- latest coord bench infra/phase1 \
  || say "warning: uncommitted changes in latest/, coord/, bench/ or infra/phase1/ are not uploaded; the instance builds $COMMIT"
L=$D/runs/$RUN
[ "$DRY" = 1 ] && L=$(mktemp -d)   # a dry run leaves no run directory
mkdir -p "$L"
TAR=$L/qap-$COMMIT.tar.gz
git -C "$R" archive --format=tar.gz -o "$TAR" HEAD Makefile latest qaplib QAPSOLN tests bench coord experiments/dnn_bound infra/phase1 \
  'experiments/estimates/*/*.par' experiments/wil50_feasibility/fullrun/fullrun.par   # the settings files recipes point at
SRC=src/qap-$COMMIT.tar.gz
say "uploading $(du -h "$TAR" | cut -f1) to s3://$B/$SRC"
[ "$DRY" = 1 ] && echo "  [dry-run] aws s3 cp $TAR s3://$B/$SRC" >&2 || aws s3 cp --only-show-errors "$TAR" "s3://$B/$SRC"
rm -f "$TAR"

# --- launch
FAM=${TYPE%%.*}
ARCH=x86_64; case $FAM in *g|*gd|*gn|*ge) ARCH=arm64 ;; esac
AMI=$(aws ssm get-parameter --name "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-$ARCH" \
  --query Parameter.Value --output text)
UD=$L/userdata_$ROLE.sh
cat > "$UD" << EOF
#!/bin/bash
export AWS_DEFAULT_REGION=$REGION QAP_RUN=$RUN QAP_BACKEND='aws:$B?region=$REGION' QAP_BUCKET=$B QAP_ROLE=$ROLE
export QAP_RECIPE=$RECIPE QAP_NO_UB=$NO_UB QAP_REFERENCE=$REFERENCE QAP_COMMIT=$COMMIT QAP_VISIBILITY=$VIS${GATES:+ QAP_GATES='$GATES'}
${NODE_ENV:+export $NODE_ENV}
shutdown -h +$LIMIT
mkdir -p /opt/qap && cd /opt/qap
for i in 1 2 3 4 5 6; do aws s3 cp --only-show-errors s3://$B/$SRC src.tgz && break; sleep 10; done
tar xzf src.tgz && bash infra/phase1/node.sh
shutdown -h now
EOF
MO=()
[ "$MARKET" = spot ] && MO=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
ID=""
for s in $SUBNETS; do   # the next zone if one has no capacity
  ID=$(change ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" --count 1 --subnet-id "$s" \
    --iam-instance-profile Name=$ROLE_NAME --security-group-ids "$SG" --instance-initiated-shutdown-behavior terminate \
    --metadata-options HttpTokens=required,HttpEndpoint=enabled \
    --block-device-mappings 'DeviceName=/dev/xvda,Ebs={VolumeSize=32,VolumeType=gp3,DeleteOnTermination=true}' \
    ${MO[@]+"${MO[@]}"} --user-data "file://$UD" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=qap-phase1-$RUN-$ROLE},{Key=qap-phase1,Value=$RUN},{Key=qap-role,Value=$ROLE}]" \
    --query 'Instances[0].InstanceId' --output text 2> "$L/launch.err") && break
  say "$TYPE: not launched in $s ($(tail -1 "$L/launch.err")); trying the next zone"
  ID=""
done
[ -n "$ID" ] || { say "$TYPE: could not launch (see $L/launch.err)"; exit 1; }
AZ=?; PRICE=?
if [ "$DRY" != 1 ]; then
  AZ=$(aws ec2 describe-instances --instance-ids "$ID" --query 'Reservations[0].Instances[0].Placement.AvailabilityZone' --output text)
  [ "$MARKET" = spot ] && PRICE=$(aws ec2 describe-spot-price-history --instance-types "$TYPE" --availability-zone "$AZ" \
    --product-descriptions Linux/UNIX --start-time "$(date -u +%FT%TZ)" --query 'SpotPriceHistory[0].SpotPrice' --output text)
fi
[ -f "$L/launches.tsv" ] || echo -e "launched\trole\ttype\tmarket\tusd_per_hour\tinstance\taz\tcommit\trecipe" > "$L/launches.tsv"
echo -e "$(date -u +%FT%TZ)\t$ROLE\t$TYPE\t$MARKET\t$PRICE\t$ID\t$AZ\t$COMMIT\t${RECIPE:--}$([ "$NO_UB" = 1 ] && echo ' no-ub')" >> "$L/launches.tsv"
echo "$REGION" > "$L/region.txt"; echo "$B" > "$L/bucket.txt"
say "$TYPE: $ID in $AZ ($MARKET, \$$PRICE/h), shuts down by itself after $LIMIT minutes"
say "build and gates take about 10-15 minutes. Then: infra/phase1/collect.sh $RUN"
