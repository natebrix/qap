#!/bin/bash
# Start a run of a QAPLIB problem on EC2, from the laptop, in one command
# (README.md, Starting a run):
#
#   infra/phase1/start.sh PROBLEM [--dry-run] [--profile P] [--region eu-north-1] [--yes]
#
# PROBLEM is a QAPLIB name: tho40, sko42, sko49, wil50, tai50b, ... The script
#  1. checks the AWS session: an IAM Identity Center (SSO) sign-in, and offers
#     `aws sso login` when there is none or it has expired. It stores no credentials:
#     every call goes through the AWS CLI and boto3 on that session (the CLI keeps its
#     own short-lived token in ~/.aws/sso/cache), and the instances use their own role;
#  2. works out the defaults and asks you to keep each one (Enter) or change it:
#     - the run's name: PROBLEM-YYMMDD (the UTC date), the name the website shows;
#     - the params file: the one recipes/defaults.txt names for PROBLEM (its next run's
#       settings), else recipes/PROBLEM.par, else the estimate's
#       (experiments/estimates/PROBLEM/PROBLEM.par);
#     - the incumbent: the best known (QAPSOLN/PROBLEM.sln) + the step of the objective's
#       grid (objstep.py, the solver's own test: 2 where every value is even, else 1): the
#       next value, so the run must find the best known itself and prove nothing is better;
#     - the instances: a coordinator on demand (c7a.2xlarge: it bounds the root and merges, and is
#       never reclaimed) and a worker group, EC2 Auto Scaling, of N workers (with the vCPU quota
#       left) of the type (c7a.16xlarge, or the problem's in recipes/defaults.txt) and its siblings
#       on the same CPU (m7a, r7a), Spot or on-demand, slots per worker (its cores, or the
#       problem's slots= where a process needs more memory), and the time limit (72 hours);
#     - if you ask for them: the data file, the solver flags (the grid's symmetry
#       file), the disk, the unit sizing and the public status;
#  3. shows the whole launch and asks before making it;
#  4. writes the run's recipe to infra/phase1/runs/RUN/recipe.sh and calls launch.sh once:
#     the coordinator, and the worker group (--workers N), which the coordinator starts once
#     the run exists and stops when it ends; the run is sized for every slot (--run-slots).
# For a RUN that exists: with no coordinator, one is launched to resume it (its group updated);
# with one, its worker group grows by N (a run from before worker groups: N --worker-only).
#
#   infra/phase1/start.sh --spec infra/phase1/specs/NAME.spec [--dry-run] [--yes]
#
#   --spec F    every answer from a committed spec file (infra/phase1/specs/README.md): no
#               questions but the last ("Launch?", which --yes answers). ARMS in it makes an
#               A/B: one run per arm, RUN-ARM, identical but for the params file
#   --no-commit leave the run's record (runs/RUN/: recipe.sh, launches.tsv, spec.txt)
#               uncommitted; by default it is committed and pushed once launched
#   --dry-run   launch.sh --dry-run for the first instance: the AWS calls that change
#               something are printed, none is made
#   --profile P the AWS CLI profile (AWS_PROFILE) to sign in with
#   --yes       take every default without asking (with --dry-run: a check)
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../.." && pwd)
P=""; DRY=0; YES=0; REGION=eu-north-1; SPEC=""; COMMIT_RECORD=1
while [ $# -gt 0 ]; do
  case $1 in
    --dry-run) DRY=1; shift ;;
    --spec) SPEC=$2; shift 2 ;;
    --no-commit) COMMIT_RECORD=0; shift ;;
    --yes) YES=1; shift ;;
    --profile) export AWS_PROFILE=$2; shift 2 ;;
    --region) REGION=$2; shift 2 ;;
    -h|--help) awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 0 ;;
    -*) echo "unknown argument $1"; exit 2 ;;
    *) [ -z "$P" ] || { echo "one PROBLEM only ($P, $1)"; exit 2; }; P=$1; shift ;;
  esac
done
[ -n "$P$SPEC" ] || { awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 2; }
[ -z "$SPEC" ] || SPEC=$(cd "$(dirname "$SPEC")" && pwd)/$(basename "$SPEC")
cd "$R"

# --- a spec file (--spec): its keys answer the questions; a key it does not set takes the default.
# It is sourced (a committed file of this repository, as recipes/*.sh), and must be committed as it
# is, so that the run's record names a file anyone can read at that commit
SPEC_KEYS="PROBLEM RUN PARAMS ARMS INCUMBENT DATA FLAGS UNITS VISIBILITY ESTIMATE ESTIMATE_ENV WORKERS WORKER_TYPE WORKER_TYPES MARKET WORKER_SLOTS COORD_TYPE COORD_SLOTS HOURS DISK"
NOASK=$YES
SPECREL=""
if [ -n "$SPEC" ]; then
  [ -f "$SPEC" ] || { echo "no spec file $SPEC"; exit 2; }
  SPECREL=${SPEC#"$R"/}
  git ls-files --error-unmatch "$SPECREL" > /dev/null 2>&1 && git diff --quiet HEAD -- "$SPECREL" \
    || { echo "$SPECREL is not committed as it is: commit it first (the run's record names it)"; exit 2; }
  BAD=$(grep -oE '^[A-Za-z_][A-Za-z0-9_]*=' "$SPEC" | tr -d = | while read -r k; do
          case " $SPEC_KEYS " in *" $k "*) ;; *) echo "$k" ;; esac; done)
  [ -z "$BAD" ] || { echo "$SPECREL: unknown key(s) $(echo $BAD) (known: $SPEC_KEYS)"; exit 2; }
  eval "$( (set +u; source "$SPEC" > /dev/null
            for k in $SPEC_KEYS; do [ -n "${!k+x}" ] && printf 'S_%s=%q\n' "$k" "${!k}"; done) )"
  [ -n "${S_PROBLEM:-}" ] || { echo "$SPECREL: no PROBLEM"; exit 2; }
  [ -z "$P" ] || [ "$P" = "$S_PROBLEM" ] || { echo "the spec is for $S_PROBLEM, not $P"; exit 2; }
  P=$S_PROBLEM; NOASK=1
  [ -z "${S_ARMS:-}" ] || [ -z "${S_PARAMS:-}" ] || { echo "$SPECREL: PARAMS or ARMS, not both"; exit 2; }
  echo "spec: $SPECREL"
fi
specval() {   # specval VAR: the spec's value for start.sh's VAR, if the spec sets it
  local k=""
  case $1 in
    RUN) k=RUN ;; PAR) k=PARAMS ;; UB) k=INCUMBENT ;; DAT) k=DATA ;; EXTRA) k=FLAGS ;; COORD) k=UNITS ;;
    VIS) k=VISIBILITY ;; ESTIMATE) k=ESTIMATE ;; N) k=WORKERS ;; TYPE) k=WORKER_TYPE ;; WTYPES) k=WORKER_TYPES ;;
    MARKET) k=MARKET ;; SLOTS) k=WORKER_SLOTS ;; CTYPE) k=COORD_TYPE ;; CSLOTS) k=COORD_SLOTS ;; HOURS) k=HOURS ;;
    DISK) k=DISK ;; P|P2) k=PROBLEM ;;
  esac
  [ -n "$k" ] && [ -n "$SPEC" ] && eval "[ -n \"\${S_$k+x}\" ]" && eval "printf '%s' \"\$S_$k\""
}

ask() {   # ask VAR QUESTION DEFAULT: Enter keeps the default; a spec's value answers it
  local v="" sv
  if sv=$(specval "$1"); then echo "$2: $sv (spec)"; v=$sv
  elif [ "$NOASK" = 1 ]; then echo "$2: $3"
  else read -r -p "$2 [$3]: " v; fi
  printf -v "$1" '%s' "${v:-$3}"
}
yesno() {   # yesno QUESTION y|n (the default): status 0 for yes
  local v=""
  if [ "$NOASK" = 1 ]; then echo "$1 $2"; v=$2; else read -r -p "$1 [$([ "$2" = y ] && echo Y/n || echo y/N)]: " v; fi
  case ${v:-$2} in y|Y|yes|Yes|YES) return 0 ;; *) return 1 ;; esac
}
mask() { if [ -n "${ACCT:-}" ]; then sed "s/$ACCT/<account>/g"; else cat; fi; }
note() { echo "  $*" | mask; }
grid() {   # STEP and GRIDWHAT for $DAT: the step of the objective's grid, by the solver's test
  local g
  g=$(python3 "$D/objstep.py" "$DAT") || { echo "cannot read $DAT"; exit 1; }
  STEP=${g%% *}; GRIDWHAT=${g#* }
  [[ $STEP =~ ^[12]$ ]] || { echo "objstep.py: $g"; exit 1; }
}
ubnote() {   # UBNOTE: what the incumbent UB means against the best known BEST
  local d
  UBNOTE=""
  [ -n "$UB" ] && [ -n "$BEST" ] || return 0
  d=$((UB - BEST))
  if [ "$d" -le 0 ]; then UBNOTE=" (not above the best known $BEST: the run cannot find it, only prove that nothing better than $UB exists)"
  elif [ "$d" -lt "$STEP" ]; then UBNOTE=" (off the objective's grid: the solver rounds it up to the next value, $BEST + $STEP, and the run must find $BEST itself)"
  else
    local past=""   # (not $( [ ] && echo ) in the assignment: under set -e its status 1 would end the script)
    if [ "$d" -gt "$STEP" ]; then past=", past the next value $BEST + $STEP"; fi
    UBNOTE=" (the best known $BEST + $d$past: the run must find $BEST itself)"
  fi
}

# --- the tools and the AWS session (SSO; nothing stored here)
command -v aws > /dev/null || { echo "the AWS CLI is not installed (brew install awscli)"; exit 1; }
python3 -c "import boto3" 2> /dev/null || { echo "python3 has no boto3: pip3 install boto3"; exit 1; }
export AWS_DEFAULT_REGION=$REGION AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
if [ -n "${AWS_ACCESS_KEY_ID:-}" ]; then
  echo "note: AWS_ACCESS_KEY_ID is set in this shell, so the AWS CLI uses those keys before any SSO sign-in."
  echo "      unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN to use your SSO sign-in."
fi
ARN=$(aws sts get-caller-identity --query Arn --output text 2> /dev/null) || ARN=""
if [ -z "$ARN" ]; then
  echo "no AWS session$([ -n "${AWS_PROFILE:-}" ] && echo " for profile $AWS_PROFILE"), or it has expired"
  if yesno "Sign in now (aws sso login)?" y; then
    command aws sso login || { echo "aws sso login failed (once per laptop: aws configure sso)"; exit 1; }
    ARN=$(aws sts get-caller-identity --query Arn --output text 2> /dev/null) || ARN=""
  fi
  [ -n "$ARN" ] || { echo "no AWS session: aws sso login, then run this again"; exit 1; }
fi
ACCT=$(echo "$ARN" | cut -d: -f5)
WHO=$(echo "$ARN" | cut -d: -f6-)
case $WHO in
  assumed-role/AWSReservedSSO_*)
    PSET=$(echo "$WHO" | cut -d/ -f2 | sed -E 's/^AWSReservedSSO_//; s/_[0-9a-f]{16}$//')
    echo "AWS: signed in with SSO (permission set $PSET), account ending ${ACCT#????????}, region $REGION" ;;
  *)
    echo "AWS: signed in as $(echo "$WHO" | cut -d/ -f1), not with SSO: this script expects an IAM Identity Center sign-in (aws sso login)" | mask
    yesno "Go on with this identity?" n || exit 1 ;;
esac
BK="aws:qap-coord-$ACCT-$REGION?region=$REGION${QAP_AWS_ENDPOINT:+&endpoint=$QAP_AWS_ENDPOINT}"   # QAP_AWS_ENDPOINT: tests only (moto)

# --- the repository: the instances build HEAD
COMMIT=$(git rev-parse --short HEAD)
BRANCH=$(git rev-parse --abbrev-ref HEAD)
git diff --quiet HEAD -- latest coord bench infra/phase1 experiments/estimates \
  || echo "note: uncommitted changes in latest/, coord/, bench/, infra/phase1/ or experiments/estimates/ are not uploaded: the instances build $COMMIT"
[ -n "$(git branch -r --contains HEAD 2> /dev/null)" ] \
  || echo "note: $COMMIT is not on origin yet; push it (git push origin $BRANCH) so the run's commit is in the repository"

# --- the problem and its defaults
while :; do
  [[ $P =~ ^[a-z]+[0-9]+[a-z]?$ ]] && [ -f "qaplib/$P.dat" ] && break
  echo "no QAPLIB problem $P (qaplib/$P.dat)"; [ "$NOASK" = 1 ] && exit 1
  ask P "Problem (a QAPLIB name)" "tho40"
done
ask P2 "Problem" "$P"
[ "$P2" = "$P" ] || { [[ $P2 =~ ^[a-z]+[0-9]+[a-z]?$ ]] && [ -f "qaplib/$P2.dat" ] || { echo "no QAPLIB problem $P2"; exit 1; }
                      ARGS=(); [ "$DRY" = 1 ] && ARGS+=(--dry-run); [ -n "${AWS_PROFILE:-}" ] && ARGS+=(--profile "$AWS_PROFILE")
                      exec "$0" "$P2" --region "$REGION" ${ARGS[@]+"${ARGS[@]}"}; }
# the best known value: the first line's second field (some .sln files end their lines with \r)
BEST=""; [ -f "QAPSOLN/$P.sln" ] && BEST=$(awk 'NR == 1 {print $2; exit}' "QAPSOLN/$P.sln" | tr -d '\r[:space:]')
[[ $BEST =~ ^[0-9]+$ ]] || { [ -z "$BEST" ] || echo "note: QAPSOLN/$P.sln does not start with n and a value; no best known"; BEST=""; }
# the data file and its symmetry: the grid-ordered copy (qaplib/Pt.dat) where its symmetry file exists
DAT=qaplib/$P.dat; EXTRA=""
if [ -f "qaplib/${P}t.dat" ] && [ -f "latest/SYM/${P}t.sym" ]; then DAT=qaplib/${P}t.dat; EXTRA="-s SYM/${P}t.sym"
elif [ -f "latest/SYM/$P.sym" ]; then EXTRA="-s SYM/$P.sym"; fi
# the params file: recipes/defaults.txt's for P (its next run's settings), else the first that exists
PAR=$(awk -v p="$P" '!/^#/ && $1 == p {print $2; exit}' infra/phase1/recipes/defaults.txt 2> /dev/null || true)
# and its instance type and slots per instance, if the line gives them (type=T slots=S after the file)
DEF_TYPE=$(awk -v p="$P" '!/^#/ && $1 == p {for (i = 3; i <= NF; i++) if ($i ~ /^type=/) {sub(/^type=/, "", $i); print $i}; exit}' infra/phase1/recipes/defaults.txt 2> /dev/null || true)
DEF_SLOTS=$(awk -v p="$P" '!/^#/ && $1 == p {for (i = 3; i <= NF; i++) if ($i ~ /^slots=/) {sub(/^slots=/, "", $i); print $i}; exit}' infra/phase1/recipes/defaults.txt 2> /dev/null || true)
# mem=G: GiB a solver process needs (tai60b: 3); the slots then default to what the type's memory holds
DEF_MEM=$(awk -v p="$P" '!/^#/ && $1 == p {for (i = 3; i <= NF; i++) if ($i ~ /^mem=/) {sub(/^mem=/, "", $i); print $i}; exit}' infra/phase1/recipes/defaults.txt 2> /dev/null || true)
[[ ${DEF_MEM:-1} =~ ^[1-9][0-9]*$ ]] || { echo "recipes/defaults.txt: mem=$DEF_MEM for $P is not whole GiB"; exit 1; }
[ -z "$PAR" ] || [ -f "$PAR" ] || { echo "note: recipes/defaults.txt names $PAR for $P, which does not exist"; PAR=""; }
if [ -z "$PAR" ]; then
  for c in "infra/phase1/recipes/$P.par" "experiments/estimates/$P/$P.par" "experiments/${P}_feasibility/fullrun/fullrun.par"; do
    [ -f "$c" ] && { PAR=$c; break; }
  done
fi
ENVV=""
COORD="--budget 1200 --target-unit-cpu 1200 --steal-after 300 --checkpoint-every 60 --ckpt-upload 300 --trace"
if [ -f "$D/recipes/$P.sh" ]; then   # a committed recipe for P: its data file, flags, environment and units
  eval "$(INSTANCE= PARAMS= EXTRA= ENV= COORD=; source "$D/recipes/$P.sh" > /dev/null
          printf 'DAT=%q EXTRA=%q ENVV=%q COORD=%q\n' "$INSTANCE" "$EXTRA" "$ENV" "$COORD")"
fi

# --- the run's name; an existing run takes instances with its own settings
MODE=new
RUN=$P-$(date -u +%y%m%d)
while :; do
  ask RUN "Run name (the website shows it)" "$RUN"
  if ! [[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]]; then echo "letters, digits, - and _ only"; [ "$NOASK" = 1 ] && exit 1; continue; fi
  ST=$(python3 "$D/bk.py" "$BK" status "$RUN" 2> /dev/null || echo none)
  case $ST in
    none) break ;;
    finished|stopped) echo "run $RUN is $ST: give the new run another name (e.g. ${RUN}b)"; [ "$NOASK" = 1 ] && exit 1; RUN=${RUN}b ;;
    *) COORD_LIVE=$(aws ec2 describe-instances --filters "Name=tag:qap-phase1,Values=$RUN" Name=tag:qap-role,Values=all \
         Name=instance-state-name,Values=pending,running --query 'Reservations[].Instances[].InstanceId' --output text) || COORD_LIVE=""
       echo "run $RUN exists ($ST), $([ -n "$COORD_LIVE" ] && echo "coordinated by $COORD_LIVE" || echo "with no instance coordinating it")"
       [ -z "$SPEC" ] || { echo "a spec starts new runs: give it another RUN (start.sh $P without --spec adds to a run)"; exit 1; }
       if yesno "Add instances to it, with its own settings?" y; then
         MODE=add; [ -n "$COORD_LIVE" ] || MODE=resume; break
       fi
       [ "$NOASK" = 1 ] && exit 1; RUN=${RUN}b ;;
  esac
done

# --- an A/B (a spec's ARMS="NAME=PARAMS NAME=PARAMS ..."): one run per arm, RUN-NAME, identical but
# for the params file (docs/review_nov18_runs.md, item 2: an A/B is not typed at prompts)
ARM_NAMES=""; ARM_PARS=""
if [ -n "${S_ARMS:-}" ]; then
  for a in $S_ARMS; do
    an=${a%%=*}; ap=${a#*=}
    [[ $an =~ ^[A-Za-z0-9_]{1,12}$ ]] && [ "$an" != "$a" ] || { echo "ARMS: $a is not NAME=PARAMS (NAME: letters, digits, _)"; exit 1; }
    [ -f "$ap" ] && git ls-files --error-unmatch "$ap" > /dev/null 2>&1 && git diff --quiet HEAD -- "$ap" \
      || { echo "ARMS: $ap is not a committed file as it is"; exit 1; }
    st=$(python3 "$D/bk.py" "$BK" status "$RUN-$an" 2> /dev/null || echo none)
    [ "$st" = none ] || { echo "run $RUN-$an exists ($st): give the spec another RUN"; exit 1; }
    ARM_NAMES="$ARM_NAMES $an"; ARM_PARS="$ARM_PARS $ap"
  done
  set -- $ARM_NAMES; [ $# -ge 2 ] || { echo "ARMS: two or more arms"; exit 1; }
  PAR=$(echo $ARM_PARS | cut -d' ' -f1)   # the first arm's, for the checks below (the others are alike)
  S_PARAMS=$PAR
fi

# --- what to solve (a new run only)
UB=""; UBNOTE=""
if [ "$MODE" = new ]; then
  while :; do
    ask PAR "Params file (relative to the repository)" "${PAR:-latest/PARAM/default.par}"
    if [ ! -f "$PAR" ]; then echo "no file $PAR"
    elif ! git ls-files --error-unmatch "$PAR" > /dev/null 2>&1 || ! git diff --quiet HEAD -- "$PAR"; then
      echo "$PAR is not committed as it is: the instances build $COMMIT, so commit it first"
    else break; fi
    [ "$NOASK" = 1 ] && exit 1
  done
  KEYS=$(grep -E '^(ship\.bits|dnn\.depths|dnn\.giveup|pool\.on|cut\.grid|cut\.zone|cut\.regions) *=' "$PAR" | sed 's/  */ /g' | paste -sd, - | sed 's/,/, /g') || KEYS=""
  note "${KEYS:-no ship.bits, dnn.depths, pool.on or cut keys: the defaults of the solver}"
  grep -qE '^ship\.bits *= *32' "$PAR" || note "ship.bits is not 32: parent states in the store take twice the room (production runs use 32)"
  grid
  note "objective values: $GRIDWHAT"
  UBDEF=none; [ -n "$BEST" ] && UBDEF=$((BEST + STEP))
  while :; do
    UBASKED=1
    ask UB "Incumbent ($([ -n "$BEST" ] && echo "the best known $BEST + $STEP, the next value; " || echo "no best known in QAPSOLN/; ")none: no incumbent)" "$UBDEF"
    [ "$UB" = none ] && break
    [[ $UB =~ ^[0-9]+$ ]] && break
    echo "a whole number, or none"; [ "$NOASK" = 1 ] && exit 1
  done
  [ "$UB" = none ] && UB=""
  ubnote
fi

# --- the instances: a coordinator on demand (the root bound and the merges: never reclaimed, so a
# run is never lost with its root; docs/review_nov18_runs.md) and a group of workers that EC2 Auto
# Scaling keeps at N, Spot by default, replacing a reclaimed one (launch.sh --workers)
TYPE=${DEF_TYPE:-c7a.16xlarge}; MARKET=spot
while :; do
  ask TYPE "Worker instance type" "$TYPE"
  CPU=$(aws ec2 describe-instance-types --instance-types "$TYPE" --query 'InstanceTypes[0].[VCpuInfo.DefaultVCpus,VCpuInfo.DefaultCores,MemoryInfo.SizeInMiB]' --output text 2> /dev/null) || CPU=""
  VCPU=$(echo "$CPU" | awk '{print $1}'); CORES=$(echo "$CPU" | awk '{print $2}'); MEM=$(echo "$CPU" | awk '{print $3}')
  [[ $MEM =~ ^[0-9]+$ ]] || MEM=
  [[ $VCPU =~ ^[1-9][0-9]*$ ]] && [[ $CORES =~ ^[1-9][0-9]*$ ]] && break
  echo "no instance type $TYPE in $REGION"; [ "$NOASK" = 1 ] && exit 1
done
# the group's types: a 7th-generation AMD type and its siblings, one CPU (AMD EPYC 9R14) in three
# families: three times the Spot pools, the same per-core speed and rounding. The slots below are
# sized by TYPE, the first, which has the least memory of the three
WTYPES=$TYPE
case $TYPE in c7a.*|m7a.*|r7a.*) WTYPES="$TYPE"; for f in c7a m7a r7a; do [ "$f" = "${TYPE%%.*}" ] || WTYPES="$WTYPES,$f.${TYPE#*.}"; done ;; esac
while :; do
  ask WTYPES "Worker group's types (comma-separated; Spot takes the cheapest with capacity)" "$WTYPES"
  [[ $WTYPES =~ ^[a-z0-9]+\.[a-z0-9]+(,[a-z0-9]+\.[a-z0-9]+)*$ ]] && break
  echo "instance types separated by commas"; [ "$NOASK" = 1 ] && exit 1
done
while :; do
  ask MARKET "Workers: spot or on-demand" "$MARKET"
  case $MARKET in spot|on-demand) break ;; esac
  echo "spot or on-demand"; [ "$NOASK" = 1 ] && exit 1
done
# the vCPU quota left: Spot "All Standard Spot Instance Requests", on-demand "Running On-Demand Standard"
quota_room() {   # quota_room MARKET VCPUS_EACH: QROOM (instances that fit), QNOTE
  local code=L-34B43A08 q used
  [ "$1" = on-demand ] && code=L-1216C47A
  q=$(aws service-quotas get-service-quota --service-code ec2 --quota-code $code --query 'Quota.Value' --output text 2> /dev/null || echo "?"); q=${q%.0}
  used=$(aws ec2 describe-instances --filters Name=instance-state-name,Values=pending,running \
    --query 'Reservations[].Instances[].[InstanceLifecycle,CpuOptions.CoreCount,CpuOptions.ThreadsPerCore]' --output text 2> /dev/null \
    | awk -v m="$1" '($1 == "spot") == (m == "spot") {s += $2 * $3} END {print s + 0}') || used=0
  [[ $used =~ ^[0-9]+$ ]] || used=0
  QROOM=""; if [[ $q =~ ^[0-9]+(\.[0-9]*)?$ ]]; then QROOM=$(( (${q%.*} - used) / $2 )); [ "$QROOM" -ge 0 ] || QROOM=0; fi
  QNOTE="$1 quota $q vCPUs, $used in use$([ -n "$QROOM" ] && echo ": room for $QROOM more of $2 vCPUs")"
}
quota_room "$MARKET" "$VCPU"; ROOM=$QROOM
note "$TYPE: $VCPU vCPUs, $CORES cores; $QNOTE"
PRICE=""
if [ "$MARKET" = spot ]; then
  PRICE=$(aws ec2 describe-spot-price-history --instance-types "$TYPE" --product-descriptions Linux/UNIX --start-time "$(date -u +%FT%TZ)" \
    --query 'SpotPriceHistory[].SpotPrice' --output text 2> /dev/null | tr -s '[:space:]' '\n' | sort -n | head -1) || PRICE=""
  [[ $PRICE =~ ^[0-9]+(\.[0-9]+)?$ ]] || PRICE=""
  [ -n "$PRICE" ] && note "Spot price now: \$$PRICE an hour for $TYPE (the cheapest zone)"
fi
NDEF=1; [ -n "$ROOM" ] && [ "$ROOM" -lt 1 ] && NDEF=0
while :; do
  ask N "Worker instances (the group keeps this many; a reclaimed one is replaced)" "$NDEF"
  if ! [[ $N =~ ^[1-9][0-9]*$ ]]; then echo "a whole number, 1 or more"; [ "$NOASK" = 1 ] && exit 1; continue; fi
  [ -z "$ROOM" ] || [ "$N" -le "$ROOM" ] && break
  echo "the quota has room for $ROOM: the group would launch only those"
  [ -n "$SPEC" ] && { note "the spec's $N: the group launches what the quota allows, and the rest when it has room"; break; }
  yesno "Ask for $N anyway?" n && break
  [ "$NOASK" = 1 ] && exit 1
done
# slots per instance: a worker process each; the cores, unless the problem needs more memory per process
# than the type has per core (recipes/defaults.txt mem=G: as many as (memory - 4 GiB) / G holds; slots=S
# fixes the number). tai60b: about 2.3 GB a process measured, 3 with the parents' states it keeps
fit_slots() {   # fit_slots CORES MEM_MIB DEFAULT: the slots the memory holds
  local s=${3:-$1} fit
  if [ -z "$DEF_SLOTS" ] && [ -n "$DEF_MEM" ] && [ -n "$2" ]; then
    fit=$(( ($2 / 1024 - 4) / DEF_MEM )); [ "$fit" -ge 1 ] || fit=1
    [ "$fit" -lt "$s" ] && s=$fit
  fi
  [ "$s" -le "$1" ] || s=$1
  echo "$s"
}
SDEF=$(fit_slots "$CORES" "$MEM" "${DEF_SLOTS:-$CORES}")
while :; do
  ask SLOTS "Slots per worker (one solver process each$([ -n "$MEM" ] && echo "; $((MEM / 1024)) GiB of memory$([ -n "$DEF_MEM" ] && echo ", $P needs about $DEF_MEM GiB a process")"))" "$SDEF"
  [[ $SLOTS =~ ^[1-9][0-9]*$ ]] && [ "$SLOTS" -le "$CORES" ] && break
  echo "a whole number from 1 to $CORES"; [ "$NOASK" = 1 ] && exit 1
done
if [ -n "$DEF_MEM" ] && [ -n "$MEM" ] && [ $((SLOTS * DEF_MEM + 4)) -gt $((MEM / 1024)) ]; then
  note "$SLOTS slots x $DEF_MEM GiB is more than the $((MEM / 1024)) GiB of $TYPE: the solver processes may run out of memory"
fi
# the coordinator: same CPU family, small, on demand; its cores but one run units (one is left to the
# coordinator's merges and uploads)
CTYPE=c7a.2xlarge; case $TYPE in c7a.*|m7a.*|r7a.*) ;; *) CTYPE=${TYPE%%.*}.2xlarge ;; esac
while :; do
  ask CTYPE "Coordinator instance type (on demand: bounds the root, merges, runs units)" "$CTYPE"
  CPU=$(aws ec2 describe-instance-types --instance-types "$CTYPE" --query 'InstanceTypes[0].[VCpuInfo.DefaultVCpus,VCpuInfo.DefaultCores,MemoryInfo.SizeInMiB]' --output text 2> /dev/null) || CPU=""
  CVCPU=$(echo "$CPU" | awk '{print $1}'); CCORES=$(echo "$CPU" | awk '{print $2}'); CMEM=$(echo "$CPU" | awk '{print $3}')
  [[ $CMEM =~ ^[0-9]+$ ]] || CMEM=
  [[ $CVCPU =~ ^[1-9][0-9]*$ ]] && [[ $CCORES =~ ^[1-9][0-9]*$ ]] && break
  echo "no instance type $CTYPE in $REGION"; [ "$NOASK" = 1 ] && exit 1
done
quota_room on-demand "$CVCPU"
note "$CTYPE: $CVCPU vCPUs, $CCORES cores; $QNOTE"
[ -z "$QROOM" ] || [ "$QROOM" -ge 1 ] || note "the on-demand quota has no room for $CTYPE: the launch will be refused"
CDEF=$(fit_slots "$CCORES" "$CMEM" "$([ "$CCORES" -gt 1 ] && echo $((CCORES - 1)) || echo 1)")
while :; do
  ask CSLOTS "Coordinator slots (one solver process each)" "$CDEF"
  [[ $CSLOTS =~ ^[0-9]+$ ]] && [ "$CSLOTS" -ge 1 ] && [ "$CSLOTS" -le "$CCORES" ] && break
  echo "a whole number from 1 to $CCORES"; [ "$NOASK" = 1 ] && exit 1
done
TOTAL=$((CSLOTS + N * SLOTS))
while :; do
  ask HOURS "Time limit in hours (the coordinator shuts down and the group goes to 0 then; a relaunch resumes the run)" 72
  [[ $HOURS =~ ^[1-9][0-9]*$ ]] && break
  echo "whole hours"; [ "$NOASK" = 1 ] && exit 1
done

# --- the rest, if asked
DISK=100; VIS=public
OTHER=n; for k in DATA FLAGS UNITS VISIBILITY ESTIMATE DISK; do eval "[ -n \"\${S_$k+x}\" ]" && OTHER=y; done
if yesno "Change the other settings (data file, solver flags, disk, units, public status, estimate)?" $OTHER; then
  if [ "$MODE" = new ]; then
    DAT0=$DAT; STEP0=$STEP
    while :; do ask DAT "Data file" "$DAT"; [ -f "$DAT" ] && break; echo "no file $DAT"; [ "$NOASK" = 1 ] && exit 1; done
    if [ "$DAT" != "$DAT0" ]; then
      grid
      if [ "$STEP" != "$STEP0" ] && [ -n "$BEST" ] && [ "$UB" = "$((BEST + STEP0))" ]; then
        UB=$((BEST + STEP)); note "$DAT: objective values $GRIDWHAT; the incumbent is now $UB"
      fi
      ubnote
    fi
    ask EXTRA "Solver flags (paths relative to latest/; none for no flags)" "${EXTRA:-none}"; [ "$EXTRA" = none ] && EXTRA=""
    ask COORD "Unit sizing (qapctl init)" "$COORD"
    ask VIS "Public or private status" "$VIS"
    # an estimate run (docs/estimates_plan.md): the solve for SOLVE_HOURS of its wall time, then
    # DIVE_SECONDS of survivor dives per open unit; the solve stays paused, resumable
    while :; do
      ask ESTIMATE "Estimate run: SOLVE_HOURS:DIVE_SECONDS[:B[:SAMPLE[:unif|qpb]]] (e.g. 2:1800, or 1:900:k1.5:0.2:qpb for a bushy tree; none: a solve)" "${ESTIMATE:-none}"
      [ "$ESTIMATE" = none ] && ESTIMATE="" && break
      [[ $ESTIMATE =~ ^[0-9]+(\.[0-9]+)?:[0-9]+(:k?[0-9.]+(:[0-9.]+(:(unif|qpb))?)?)?$ ]] && break
      echo "SOLVE_HOURS:DIVE_SECONDS[:B[:SAMPLE[:unif|qpb]]], or none"; [ "$NOASK" = 1 ] && exit 1
    done
  fi
  ask DISK "Disk per instance, GB" "$DISK"
fi
EST_ENV=""
if [ -n "${ESTIMATE:-}" ]; then
  IFS=: read -r EH ED EREST <<< "$ESTIMATE"
  EST_ENV="QAP_ESTIMATE=$(python3 -c "print(int(float('$EH') * 3600))"):$ED${EREST:+:$EREST}"
fi
# a spec's ESTIMATE_ENV: more estimator settings, QAP_EST_X=V separated by spaces (no space in a value)
EST_MORE=${S_ESTIMATE_ENV:-}
if [ -n "$EST_MORE" ]; then
  [ -n "$EST_ENV" ] || { echo "ESTIMATE_ENV without ESTIMATE"; exit 1; }
  for kv in $EST_MORE; do [[ $kv =~ ^QAP_EST_[A-Z_]+=[^[:space:]]+$ ]] || { echo "ESTIMATE_ENV: $kv is not QAP_EST_X=V"; exit 1; }; done
fi

# --- the launch, shown in full (a new run: the params file and, with a best known, an answer for the incumbent)
# (bash skips the rest of a block after an expansion error: Nov 16, a \r in tai50b.sln left it with none)
[ "$MODE" != new ] || { [ -n "$PAR" ] && [ "${UBASKED:-0}" = 1 ]; } \
  || { echo "the params file or the incumbent was not settled: nothing launched"; exit 1; }
echo
echo "Run $RUN: $([ "$MODE" = new ] && echo "a new run" || echo "add instances to the run ($MODE)"), commit $COMMIT ($BRANCH)"
if [ "$MODE" = new ]; then
  note "problem    $P: $DAT${EXTRA:+, $EXTRA}; objective values $GRIDWHAT"
  if [ -n "$ARM_NAMES" ]; then
    i=1; for an in $ARM_NAMES; do note "arm        $RUN-$an: $(echo $ARM_PARS | cut -d' ' -f$i)"; i=$((i + 1)); done
    note "           the arms differ only in their params files:"
    # (diff's status is 1 when the files differ: under pipefail the pipeline must not end the script)
    set -- $ARM_PARS; for ap in "${@:2}"; do { diff -u "$1" "$ap" || true; } | sed -n '3,$p' | grep -E '^[-+]' | sed 's/^/             /' | head -20 || true; done
  else
    note "params     $PAR"
  fi
  note "incumbent  ${UB:-none}$UBNOTE"
  note "units      $COORD"
  note "status     $VIS"
  [ -n "$EST_ENV" ] && note "estimate   the solve for ${EH} h of its wall time, then $ED s of dives a unit as run $RUN-dives (${EST_ENV}); the solve stays paused"
  [ -n "$EST_MORE" ] && note "           and $EST_MORE"
fi
note "coordinator $CTYPE, on demand, $REGION: $CSLOTS slots ($CCORES cores$([ -n "$CMEM" ] && echo ", $((CMEM / 1024)) GiB"))"
note "workers    a group of $N x ($WTYPES), $MARKET, every zone of $REGION: $SLOTS slots each ($CORES cores$([ -n "$MEM" ] && echo ", $((MEM / 1024)) GiB"))$([ "$MODE" = new ] && echo "; the run sized for $TOTAL slots")"
note "           started once the run exists (the root bounded), replaced when Spot reclaims one, 0 when the run ends"
[ -n "$PRICE" ] && note "cost       workers about \$$(python3 -c "print(f'{$PRICE * $N * 24:.0f}')") a day at today's Spot price, plus the coordinator on demand"
note "limit      $HOURS hours ($((HOURS * 60)) minutes); disk $DISK GB each"
[ "$DRY" = 1 ] && note "dry run: launch.sh --dry-run for the first instance; nothing is created"
echo
[ -z "$SPEC" ] || note "the spec    $SPECREL$([ -n "$ARM_NAMES" ] && echo ", $(echo $ARM_NAMES | wc -w | tr -d ' ') runs")"
[ "$COMMIT_RECORD" = 1 ] && [ "$DRY" != 1 ] && [ "$MODE" = new ] && note "record     runs/$RUN$([ -n "$ARM_NAMES" ] && echo "-*")/ committed and pushed once launched (--no-commit: not)"
echo
# the last question is asked even with a spec; --yes answers it
if [ "$YES" = 1 ]; then echo "Launch? y"
else read -r -p "Launch? [y/N]: " v; case $v in y|Y|yes|Yes|YES) ;; *) echo "nothing launched"; exit 0 ;; esac; fi

# --- the recipe (a new run) and the launch: the coordinator, and with it the worker group
COMMON=(--region "$REGION" --limit $((HOURS * 60)) --disk "$DISK" --no-reference)
[ "$DRY" = 1 ] && COMMON+=(--dry-run)
CENV="QAP_SLOTS=$CSLOTS"; [ -n "$EST_ENV" ] && CENV="$CENV $EST_ENV${EST_MORE:+ $EST_MORE}"   # node.sh: the coordinator's slots; an estimate's flow and settings
GROUP=(--workers "$N" --worker-types "$WTYPES"); [ "$SLOTS" = "$CORES" ] || GROUP+=(--worker-slots "$SLOTS")
[ "$MARKET" = on-demand ] && GROUP+=(--workers-on-demand)
[ -n "$EST_ENV" ] && GROUP+=(--worker-env "$EST_ENV")
case $MODE in
  new)
    # one run, or one per arm of an A/B (RUN-ARM, each with its own coordinator and worker group)
    RUNS=$RUN; PARS=$PAR
    if [ -n "$ARM_NAMES" ]; then RUNS=""; for an in $ARM_NAMES; do RUNS="$RUNS $RUN-$an"; done; PARS=$ARM_PARS; fi
    set -- $PARS
    for R1 in $RUNS; do
      P1=$1; shift
      RD=$D/runs/$R1; [ "$DRY" = 1 ] && RD=$(mktemp -d)   # a dry run leaves no run directory
      mkdir -p "$RD"; RF=$RD/recipe.sh
      PL=../$P1; case $P1 in latest/*) PL=${P1#latest/} ;; esac   # PARAMS is relative to latest/
      { echo "# run $R1: written by infra/phase1/start.sh at $(date -u +%FT%TZ), commit $COMMIT${SPECREL:+, spec $SPECREL}"
        echo "# $P, the best known ${BEST:-?} (QAPSOLN/$P.sln). INSTANCE is relative to the repository,"
        echo "# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh"
        printf 'INSTANCE=%q\nPARAMS=%q\nOPT=%q\nUB=%q\nEXTRA=%q\nENV=%q\nCOORD=%q\n' "$DAT" "$PL" "$BEST" "$UB" "$EXTRA" "$ENVV" "$COORD"
      } > "$RF"
      # what was launched, and from what: the spec as committed, and the answers as taken
      { [ -n "$SPECREL" ] && { echo "# spec $SPECREL at $COMMIT"; cat "$SPEC"; echo; }
        echo "# launched $(date -u +%FT%TZ): coordinator $CTYPE on demand, $CSLOTS slots; workers $N x ($WTYPES) $MARKET,"
        echo "# $SLOTS slots each; $TOTAL slots; limit $HOURS h; disk $DISK GB; status $VIS${EST_ENV:+; $EST_ENV}${EST_MORE:+ $EST_MORE}"
      } > "$RD/spec.txt"
      echo "recipe: $RF"
      FIRST=(--recipe-file "$RF" --run-slots "$TOTAL"); [ "$VIS" = private ] && FIRST+=(--private)
      echo "--- $R1: the coordinator and the worker group"
      "$D/launch.sh" "$R1" "${FIRST[@]}" --type "$CTYPE" --on-demand "${GROUP[@]}" "${COMMON[@]}" --node-env "$CENV" \
        && LAUNCHED="${LAUNCHED:-} $R1" || echo "$R1 was not launched (see above)"
    done ;;
  resume)
    echo "--- the coordinator (resuming $RUN) and the worker group"
    "$D/launch.sh" "$RUN" --type "$CTYPE" --on-demand "${GROUP[@]}" "${COMMON[@]}" --node-env "$CENV" ;;
  add)
    # a coordinator runs it: its worker group grows by N, or (a run from before worker groups) N workers
    if [ "$(aws autoscaling describe-auto-scaling-groups --auto-scaling-group-names "qap-$RUN" --query 'length(AutoScalingGroups)' --output text 2> /dev/null)" = 1 ]; then
      read -r GD GM <<< "$(aws autoscaling describe-auto-scaling-groups --auto-scaling-group-names "qap-$RUN" --query 'AutoScalingGroups[0].[DesiredCapacity,MaxSize]' --output text)"
      echo "--- the worker group qap-$RUN: $GD wanted, $((GD + N)) from now"
      if [ "$DRY" = 1 ]; then echo "  [dry-run] aws autoscaling update-auto-scaling-group --auto-scaling-group-name qap-$RUN --max-size $((GM > GD + N ? GM : GD + N)) --desired-capacity $((GD + N))"
      else aws autoscaling update-auto-scaling-group --auto-scaling-group-name "qap-$RUN" --max-size $((GM > GD + N ? GM : GD + N)) --desired-capacity $((GD + N)); fi
    else
      WENV=""; [ "$SLOTS" = "$CORES" ] || WENV="QAP_SLOTS=$SLOTS"; [ -n "$EST_ENV" ] && WENV="${WENV:+$WENV }$EST_ENV"
      W=(--type "$TYPE"); [ "$MARKET" = on-demand ] && W+=(--on-demand)
      i=1   # not seq: BSD seq counts down from 2 to 1
      while [ "$i" -le "$N" ]; do
        echo "--- worker $i of $N"
        "$D/launch.sh" "$RUN" --worker-only "${W[@]}" "${COMMON[@]}" ${WENV:+--node-env "$WENV"} || echo "worker $i was not launched"
        i=$((i + 1))
      done
    fi ;;
esac
# --- the record, committed and pushed at once (docs/review_nov18_runs.md, item 5: a run whose record
# is on one laptop only can be collected from nowhere else); only the runs' own files, which
# .gitignore keeps clear of the bucket's name (bucket.txt, userdata_*.sh)
if [ -n "${LAUNCHED:-}" ] && [ "$DRY" != 1 ]; then
  PATHS=""; for R1 in $LAUNCHED; do PATHS="$PATHS infra/phase1/runs/$R1"; done
  if [ "$COMMIT_RECORD" = 1 ]; then
    git add -- $PATHS && git commit -q -m "phase 1 runs/$(echo $LAUNCHED | sed 's/ /, /g'): launched${SPECREL:+ from $SPECREL} at $COMMIT" -- $PATHS \
      && echo "record: committed $(git rev-parse --short HEAD)" || echo "record: not committed (git add/commit failed); commit $PATHS yourself"
    git push -q origin HEAD 2> /dev/null && echo "record: pushed to origin/$BRANCH" \
      || echo "record: not pushed (git pull --no-rebase origin $BRANCH && git push)"
  else
    echo "Commit the run's record: git add$PATHS && git commit -m 'launched' && git push"
  fi
fi
echo
for R1 in ${LAUNCHED:-$RUN}; do echo "Next: infra/phase1/collect.sh $R1 (any time: instances, status, steps, the readout once it is done)."; done
echo "The coordinator's build and gates take 10-15 minutes; its workers start once the run exists."
