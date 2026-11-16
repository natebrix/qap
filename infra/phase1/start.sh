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
#     - the params file: recipes/PROBLEM.par if there is one, else the estimate's
#       (experiments/estimates/PROBLEM/PROBLEM.par);
#     - the incumbent: the best known (QAPSOLN/PROBLEM.sln) + the step of the objective's
#       grid (objstep.py, the solver's own test: 2 where every value is even, else 1): the
#       next value, so the run must find the best known itself and prove nothing is better;
#     - the instances: how many (with the vCPU quota left), the type (c7a.16xlarge),
#       Spot or on-demand, and the time limit (72 hours);
#     - if you ask for them: the data file, the solver flags (the grid's symmetry
#       file), the disk, the unit sizing and the public status;
#  3. shows the whole launch and asks before making it;
#  4. writes the run's recipe to infra/phase1/runs/RUN/recipe.sh and calls launch.sh
#     for each instance: the first coordinates the run, the others are workers
#     (--worker-only), and the run is sized for all of their cores (--run-slots).
# For a RUN that exists, the instances are added to it, with its own settings (the
# first one coordinates it if no instance does).
#
#   --dry-run   launch.sh --dry-run for the first instance: the AWS calls that change
#               something are printed, none is made
#   --profile P the AWS CLI profile (AWS_PROFILE) to sign in with
#   --yes       take every default without asking (with --dry-run: a check)
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../.." && pwd)
P=""; DRY=0; YES=0; REGION=eu-north-1
while [ $# -gt 0 ]; do
  case $1 in
    --dry-run) DRY=1; shift ;;
    --yes) YES=1; shift ;;
    --profile) export AWS_PROFILE=$2; shift 2 ;;
    --region) REGION=$2; shift 2 ;;
    -h|--help) awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 0 ;;
    -*) echo "unknown argument $1"; exit 2 ;;
    *) [ -z "$P" ] || { echo "one PROBLEM only ($P, $1)"; exit 2; }; P=$1; shift ;;
  esac
done
[ -n "$P" ] || { awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 2; }
cd "$R"

ask() {   # ask VAR QUESTION DEFAULT: Enter keeps the default
  local v=""
  if [ "$YES" = 1 ]; then echo "$2: $3"; else read -r -p "$2 [$3]: " v; fi
  printf -v "$1" '%s' "${v:-$3}"
}
yesno() {   # yesno QUESTION y|n (the default): status 0 for yes
  local v=""
  if [ "$YES" = 1 ]; then echo "$1 $2"; v=$2; else read -r -p "$1 [$([ "$2" = y ] && echo Y/n || echo y/N)]: " v; fi
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
  elif [ "$d" -lt "$STEP" ]; then UBNOTE=" (below the next value, $BEST + $STEP: the solver fathoms within its tolerance, so the run cannot find $BEST, only prove that nothing is better)"
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
  echo "no QAPLIB problem $P (qaplib/$P.dat)"; [ "$YES" = 1 ] && exit 1
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
PAR=""
for c in "infra/phase1/recipes/$P.par" "experiments/estimates/$P/$P.par" "experiments/${P}_feasibility/fullrun/fullrun.par"; do
  [ -f "$c" ] && { PAR=$c; break; }
done
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
  if ! [[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]]; then echo "letters, digits, - and _ only"; [ "$YES" = 1 ] && exit 1; continue; fi
  ST=$(python3 "$D/bk.py" "$BK" status "$RUN" 2> /dev/null || echo none)
  case $ST in
    none) break ;;
    finished|stopped) echo "run $RUN is $ST: give the new run another name (e.g. ${RUN}b)"; [ "$YES" = 1 ] && exit 1; RUN=${RUN}b ;;
    *) COORD_LIVE=$(aws ec2 describe-instances --filters "Name=tag:qap-phase1,Values=$RUN" Name=tag:qap-role,Values=all \
         Name=instance-state-name,Values=pending,running --query 'Reservations[].Instances[].InstanceId' --output text) || COORD_LIVE=""
       echo "run $RUN exists ($ST), $([ -n "$COORD_LIVE" ] && echo "coordinated by $COORD_LIVE" || echo "with no instance coordinating it")"
       if yesno "Add instances to it, with its own settings?" y; then
         MODE=add; [ -n "$COORD_LIVE" ] || MODE=resume; break
       fi
       [ "$YES" = 1 ] && exit 1; RUN=${RUN}b ;;
  esac
done

# --- what to solve (a new run only)
UB=""; UBNOTE=""
if [ "$MODE" = new ]; then
  while :; do
    ask PAR "Params file (relative to the repository)" "${PAR:-latest/PARAM/default.par}"
    if [ ! -f "$PAR" ]; then echo "no file $PAR"
    elif ! git ls-files --error-unmatch "$PAR" > /dev/null 2>&1 || ! git diff --quiet HEAD -- "$PAR"; then
      echo "$PAR is not committed as it is: the instances build $COMMIT, so commit it first"
    else break; fi
    [ "$YES" = 1 ] && exit 1
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
    echo "a whole number, or none"; [ "$YES" = 1 ] && exit 1
  done
  [ "$UB" = none ] && UB=""
  ubnote
fi

# --- the instances
TYPE=c7a.16xlarge; MARKET=spot
while :; do
  ask TYPE "Instance type" "$TYPE"
  CPU=$(aws ec2 describe-instance-types --instance-types "$TYPE" --query 'InstanceTypes[0].VCpuInfo.[DefaultVCpus,DefaultCores]' --output text 2> /dev/null) || CPU=""
  VCPU=$(echo "$CPU" | awk '{print $1}'); CORES=$(echo "$CPU" | awk '{print $2}')
  [[ $VCPU =~ ^[1-9][0-9]*$ ]] && [[ $CORES =~ ^[1-9][0-9]*$ ]] && break
  echo "no instance type $TYPE in $REGION"; [ "$YES" = 1 ] && exit 1
done
while :; do
  ask MARKET "Spot or on-demand" "$MARKET"
  case $MARKET in spot|on-demand) break ;; esac
  echo "spot or on-demand"; [ "$YES" = 1 ] && exit 1
done
# the vCPU quota left: Spot "All Standard Spot Instance Requests", on-demand "Running On-Demand Standard"
QCODE=L-34B43A08; [ "$MARKET" = on-demand ] && QCODE=L-1216C47A
Q=$(aws service-quotas get-service-quota --service-code ec2 --quota-code $QCODE --query 'Quota.Value' --output text 2> /dev/null || echo "?"); Q=${Q%.0}
USED=$(aws ec2 describe-instances --filters Name=instance-state-name,Values=pending,running \
  --query 'Reservations[].Instances[].[InstanceLifecycle,CpuOptions.CoreCount,CpuOptions.ThreadsPerCore]' --output text 2> /dev/null \
  | awk -v m="$MARKET" '($1 == "spot") == (m == "spot") {s += $2 * $3} END {print s + 0}') || USED=0
[[ $USED =~ ^[0-9]+$ ]] || USED=0
ROOM=""; if [[ $Q =~ ^[0-9]+(\.[0-9]*)?$ ]]; then ROOM=$(( (${Q%.*} - USED) / VCPU )); [ "$ROOM" -ge 0 ] || ROOM=0; fi
note "$TYPE: $VCPU vCPUs, $CORES cores; $MARKET quota $Q vCPUs, $USED in use$([ -n "$ROOM" ] && echo ": room for $ROOM more")"
PRICE=""
if [ "$MARKET" = spot ]; then
  PRICE=$(aws ec2 describe-spot-price-history --instance-types "$TYPE" --product-descriptions Linux/UNIX --start-time "$(date -u +%FT%TZ)" \
    --query 'SpotPriceHistory[].SpotPrice' --output text 2> /dev/null | tr -s '[:space:]' '\n' | sort -n | head -1) || PRICE=""
  [[ $PRICE =~ ^[0-9]+(\.[0-9]+)?$ ]] || PRICE=""
  [ -n "$PRICE" ] && note "Spot price now: \$$PRICE an hour (the cheapest zone)"
fi
NDEF=1; [ -n "$ROOM" ] && [ "$ROOM" -lt 1 ] && NDEF=0
while :; do
  ask N "Number of instances" "$NDEF"
  if ! [[ $N =~ ^[1-9][0-9]*$ ]]; then echo "a whole number, 1 or more"; [ "$YES" = 1 ] && exit 1; continue; fi
  [ -z "$ROOM" ] || [ "$N" -le "$ROOM" ] && break
  echo "the quota has room for $ROOM: the launches past that would be refused"
  yesno "Launch $N anyway?" n && break
  [ "$YES" = 1 ] && exit 1
done
TOTAL=$((N * CORES))
while :; do
  ask HOURS "Time limit in hours (each instance shuts itself down then; a relaunch resumes the run)" 72
  [[ $HOURS =~ ^[1-9][0-9]*$ ]] && break
  echo "whole hours"; [ "$YES" = 1 ] && exit 1
done

# --- the rest, if asked
DISK=100; VIS=public
if yesno "Change the other settings (data file, solver flags, disk, units, public status)?" n; then
  if [ "$MODE" = new ]; then
    DAT0=$DAT; STEP0=$STEP
    while :; do ask DAT "Data file" "$DAT"; [ -f "$DAT" ] && break; echo "no file $DAT"; done
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
  fi
  ask DISK "Disk per instance, GB" "$DISK"
fi

# --- the launch, shown in full (a new run: the params file and, with a best known, an answer for the incumbent)
# (bash skips the rest of a block after an expansion error: Nov 16, a \r in tai50b.sln left it with none)
[ "$MODE" != new ] || { [ -n "$PAR" ] && [ "${UBASKED:-0}" = 1 ]; } \
  || { echo "the params file or the incumbent was not settled: nothing launched"; exit 1; }
echo
echo "Run $RUN: $([ "$MODE" = new ] && echo "a new run" || echo "add instances to the run ($MODE)"), commit $COMMIT ($BRANCH)"
if [ "$MODE" = new ]; then
  note "problem    $P: $DAT${EXTRA:+, $EXTRA}; objective values $GRIDWHAT"
  note "params     $PAR"
  note "incumbent  ${UB:-none}$UBNOTE"
  note "units      $COORD"
  note "status     $VIS"
fi
note "instances  $N x $TYPE, $MARKET, $REGION: $CORES cores each$([ "$MODE" = new ] && echo ", the run sized for $TOTAL slots")"
[ -n "$PRICE" ] && note "cost       about \$$(python3 -c "print(f'{$PRICE * $N * 24:.0f}')") a day at today's Spot price"
note "limit      $HOURS hours ($((HOURS * 60)) minutes); disk $DISK GB each"
[ "$DRY" = 1 ] && note "dry run: launch.sh --dry-run for the first instance; nothing is created"
echo
yesno "Launch?" "$([ "$YES" = 1 ] && echo y || echo n)" || { echo "nothing launched"; exit 0; }

# --- the recipe (a new run) and the launches
COMMON=(--type "$TYPE" --region "$REGION" --limit $((HOURS * 60)) --disk "$DISK" --no-reference)
[ "$MARKET" = on-demand ] && COMMON+=(--on-demand)
[ "$DRY" = 1 ] && COMMON+=(--dry-run)
FIRST=()
case $MODE in
  new)
    RD=$D/runs/$RUN; [ "$DRY" = 1 ] && RD=$(mktemp -d)   # a dry run leaves no run directory
    mkdir -p "$RD"; RF=$RD/recipe.sh
    PL=../$PAR; case $PAR in latest/*) PL=${PAR#latest/} ;; esac   # PARAMS is relative to latest/
    { echo "# run $RUN: written by infra/phase1/start.sh at $(date -u +%FT%TZ), commit $COMMIT"
      echo "# $P, the best known ${BEST:-?} (QAPSOLN/$P.sln). INSTANCE is relative to the repository,"
      echo "# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh"
      printf 'INSTANCE=%q\nPARAMS=%q\nOPT=%q\nUB=%q\nEXTRA=%q\nENV=%q\nCOORD=%q\n' "$DAT" "$PL" "$BEST" "$UB" "$EXTRA" "$ENVV" "$COORD"
    } > "$RF"
    echo "recipe: $RF"
    FIRST=(--recipe-file "$RF" --run-slots "$TOTAL"); [ "$VIS" = private ] && FIRST+=(--private) ;;
  add) FIRST=(--worker-only) ;;
  resume) ;;
esac
echo "--- instance 1 of $N: $([ "$MODE" = add ] && echo worker || echo "coordinator and worker")"
"$D/launch.sh" "$RUN" ${FIRST[@]+"${FIRST[@]}"} "${COMMON[@]}"
i=2   # not seq: BSD seq counts down from 2 to 1
while [ "$i" -le "$N" ]; do
  if [ "$DRY" = 1 ]; then echo "--- instance $i of $N: launch.sh $RUN --worker-only ${COMMON[*]} (not run in a dry run: it needs the first)"
  else
    echo "--- instance $i of $N: worker"
    "$D/launch.sh" "$RUN" --worker-only "${COMMON[@]}" || echo "instance $i was not launched; add it later: infra/phase1/start.sh $P (run $RUN)"
  fi
  i=$((i + 1))
done
echo
echo "Next: infra/phase1/collect.sh $RUN (any time: instances, status, steps). The build and gates take 10-15 minutes;"
echo "the website lists the run once its coordinator publishes, every 5 minutes after that."
[ "$DRY" = 1 ] || echo "Commit the run's record: git add infra/phase1/runs/$RUN && git commit -m '$RUN: launched' && git push"
