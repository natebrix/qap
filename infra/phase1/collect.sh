#!/bin/bash
# Phase 1 (README.md): where a run stands, and its results once it is done.
#   infra/phase1/collect.sh [RUN]     (default: the latest run under runs/; a run launched
#                                     elsewhere is fetched into a new runs/RUN/)
# Shows the run's instances, the coordinator's status.json, each host's last
# steps, and results.txt (or results_partial.txt: a solve paused for an estimate, or stopped), compare.txt
# and an estimate run's estimate.txt when they exist. Everything it fetches
# goes to runs/RUN/ (each host's out/ under hosts/HOST/). Safe to run at any time.
# It ends with the run's readout (runs/RUN/readout.md, coord/readout.py: where the run stands and what
# it found) and a one-line verdict on the collect itself, also in runs/RUN/collect.txt: complete, or
# what is missing (docs/review_nov18_runs.md, item 4: a missing trace said only in passing).
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
MISSING=()   # what this collect could not fetch, for the verdict
# an SSO session that ends part-way leaves a collect half done (tho40-est, Nov 18: the traces step
# failed on an expired token): with under 30 minutes left, sign in again first
python3 - << 'EOF' || exit 1
import glob, json, os, sys, time
from datetime import datetime
if os.environ.get("AWS_ACCESS_KEY_ID") or os.environ.get("AWS_WEB_IDENTITY_TOKEN_FILE"):
    sys.exit(0)                      # keys or GitHub's OIDC: no SSO session to run out
left = None
for f in glob.glob(os.path.expanduser("~/.aws/sso/cache/*.json")):
    try:
        d = json.load(open(f))
        if "accessToken" in d and "expiresAt" in d:
            t = datetime.fromisoformat(d["expiresAt"].replace("Z", "+00:00").replace("UTC", "+00:00")).timestamp()
            left = max(left or 0, t - time.time())
    except (ValueError, OSError, KeyError):
        pass
if left is not None and left < 1800:
    print(f"the AWS sign-in ends in {max(left, 0) / 60:.0f} minutes, and a collect can take an hour: "
          f"aws sso login, then run this again")
    sys.exit(1)
EOF
RUN=${1:-$(ls -1t "$D/runs" 2> /dev/null | head -1)}
L=$D/runs/$RUN
[ -n "$RUN" ] || { echo "no run under $D/runs: give its name"; exit 1; }
[[ $RUN =~ ^[A-Za-z0-9_-]{1,60}$ ]] || { echo "RUN may have letters, digits, - and _ only"; exit 1; }
# a run launched elsewhere (another clone, or the GitHub workflow) has no directory here
# yet: its region is QAP_REGION's (default eu-north-1), and the store must know the run
REGION=${QAP_REGION:-eu-north-1}
[ -f "$L/region.txt" ] && REGION=$(cat "$L/region.txt")
# bucket.txt stays out of git (the bucket's name holds the account number); without it,
# the bucket is the signed-in account's, as launch.sh names it
if [ -f "$L/bucket.txt" ]; then
  B=$(cat "$L/bucket.txt")
else
  ACCT=$(command aws sts get-caller-identity --query Account --output text) || { echo "not signed in: aws sso login"; exit 1; }
  B=qap-coord-$ACCT-$REGION
fi
if [ ! -d "$L" ]; then
  ST=$(python3 "$D/bk.py" "aws:$B?region=$REGION" status "$RUN" 2> /dev/null || echo none)
  [ "$ST" != none ] || { echo "no run $RUN in s3://qap-coord-<account>-$REGION (QAP_REGION for another region)"; exit 1; }
  mkdir -p "$L"
  echo "$REGION" > "$L/region.txt"
  echo "run $RUN ($ST) was launched elsewhere: its records go to $L"
fi
export AWS_PAGER=""
aws() { command aws --region "$REGION" "$@"; }
BKS="aws:$B?region=$REGION${QAP_AWS_ENDPOINT:+&endpoint=$QAP_AWS_ENDPOINT}"   # QAP_AWS_ENDPOINT: tests only (moto)
echo "instances of $RUN:"
# with why a stopped one stopped (StateReason: a Spot reclaim reads Server.SpotInstanceTermination,
# so an interruption is not taken for a hang; docs/review_nov18_runs.md, item 3)
aws ec2 describe-instances --filters "Name=tag:qap-phase1,Values=$RUN" \
  --query 'Reservations[].Instances[].[Tags[?Key==`qap-role`]|[0].Value,InstanceId,InstanceType,InstanceLifecycle||`on-demand`,State.Name,LaunchTime,StateReason.Code||``]' \
  --output text | sed 's/^/  /'
# the worker group (launch.sh --workers): its size, and what it did lately (Spot workers it replaced)
if [ "$(aws autoscaling describe-auto-scaling-groups --auto-scaling-group-names "qap-$RUN" \
        --query 'length(AutoScalingGroups)' --output text 2> /dev/null)" = 1 ]; then
  aws autoscaling describe-auto-scaling-groups --auto-scaling-group-names "qap-$RUN" \
    --query 'AutoScalingGroups[0].[DesiredCapacity,MaxSize,length(Instances[?LifecycleState==`InService`])]' --output text \
    | awk -v g="qap-$RUN" '{print "worker group " g ": " $1 " wanted (at most " $2 "), " $3 " in service"}'
  aws autoscaling describe-scaling-activities --auto-scaling-group-name "qap-$RUN" --max-items 6 \
    --query 'Activities[].[StartTime,StatusCode,Description]' --output text 2> /dev/null | cut -c1-160 | sed 's/^/  /'
  N_SPOT=$(aws autoscaling describe-scaling-activities --auto-scaling-group-name "qap-$RUN" --max-items 200 \
    --query 'length(Activities[?contains(Cause, `Spot`) || contains(Description, `interruption`)])' --output text 2> /dev/null || echo 0)
  [ "${N_SPOT:-0}" != 0 ] && echo "  Spot interruptions the group answered: $N_SPOT"
fi
get() { aws s3 cp --only-show-errors "s3://$B/runs/$RUN/$1" "$L/$1" 2> /dev/null; }
if get status.json; then
  python3 - "$L/status.json" << 'EOF'
import json, sys, time
s = json.load(open(sys.argv[1]))
age = time.time() - s["updated"]
print(f"status: {s['status']}, {s['wall_s'] / 60:.1f} min, updated {age / 60:.1f} min ago; units {s['units']}, "
      f"queued {s['queue']}; {s['nodes']} nodes ({s['dnn_nodes']} DNN, {s['dnn_fathomed']} DNN-fathomed), "
      f"CPU {s['cpu_s']:.0f} s, incumbent {s['incumbent']}")
EOF
else
  echo "status: no status.json yet (the instance builds and checks the solver first, about 10-15 minutes)"
fi
for h in $(aws s3 ls "s3://$B/runs/$RUN/hosts/" | awk '/PRE/ {print $2}' | tr -d /); do
  mkdir -p "$L/hosts/$h"
  aws s3 cp --only-show-errors "s3://$B/runs/$RUN/hosts/$h/out.tgz" "$L/hosts/$h/out.tgz" \
    && tar xzf "$L/hosts/$h/out.tgz" -C "$L/hosts/$h" && rm -f "$L/hosts/$h/out.tgz"
  echo "host $h, last steps:"
  tail -3 "$L/hosts/$h/out/steps.txt" | cut -c1-300 | sed 's/^/  /'
done
get results.txt && { echo; cat "$L/results.txt"; }
# a solve its operation left unfinished (paused for an estimate, stopped): its RESULTS block as far as it
# went (coord/coordinator.py, write_partial)
get results_partial.txt && { echo; cat "$L/results_partial.txt"; }
# an estimate run (QAP_ESTIMATE; coord/dives.py): its estimate is its dives run's results
if aws s3 cp --only-show-errors "s3://$B/runs/$RUN-dives/results.txt" "$L/estimate.txt" 2> /dev/null; then
  aws s3 cp --only-show-errors "s3://$B/runs/$RUN-dives/estimate.json" "$L/estimate.json" 2> /dev/null
  echo; echo "estimate (run $RUN-dives; $RUN stays paused, resumable):"; cat "$L/estimate.txt"
fi
# an estimate run's dives, for analysis (dives.jsonl.gz and units.tsv under runs/RUN/dives/): part-way too;
# a dives run of its own (a re-dive of a paused frontier: start.sh FROM= DIVES=), its results and dives
DV=$RUN-dives
if [ "$(python3 "$D/bk.py" "$BKS" kind "$RUN" 2> /dev/null || echo none)" = dives ]; then
  DV=$RUN
  get estimate.txt && get estimate.json && { echo; echo "estimate (dives run $RUN):"; cat "$L/estimate.txt"; }
fi
# a sampled run (coord/sample.py): its sample, for experiments/estimates_nov19/paired.py
get sample.json && echo "sample: $(python3 -c "import json,sys; s=json.load(open(sys.argv[1])); print(len(s['nodes']), 'open nodes of', s['meta']['src'], 'by', s['meta']['by'])" "$L/sample.json")"
if [ "$(python3 "$D/bk.py" "$BKS" status "$DV" 2> /dev/null || echo none)" != none ]; then
  if (cd "$D/../.." && python3 -m coord.qapctl dives-export "$DV" --backend "$BKS" --out "$L/dives") > "$L/dives.log" 2>&1; then
    echo "dives: $(tail -1 "$L/dives.log" | sed 's/^[0-9:]* //')"
  else
    echo "dives: $(tail -1 "$L/dives.log")"; MISSING+=("dives")
  fi
  rm -f "$L/dives.log"
fi
# a --trace run's node and DNN traces (qapctl export, coord/export.py) go to QAP_EXPORT_DIR/RUN
# (default infra/phase1/exports/RUN, kept out of git; the GitHub workflow keeps them as the job's
# artifact), whatever their size. Finished or not: part-way (still running, or ended by its
# --limit, which writes no results.txt) the export holds the units merged so far, as its
# manifest says; a run without --trace has nothing to export
EXPORTS=${QAP_EXPORT_DIR:-$D/exports}
if (cd "$D/../.." && python3 -m coord.qapctl export "$RUN" --backend "$BKS" \
      --out "$EXPORTS/$RUN") > "$EXPORTS.$RUN.log" 2>&1; then
  echo "traces: $(tail -1 "$EXPORTS.$RUN.log" | sed 's/^[0-9:]* //')"
else
  if ! grep -q "nothing to export" "$EXPORTS.$RUN.log"; then
    echo "traces: $(tail -1 "$EXPORTS.$RUN.log")"; MISSING+=("traces")
  fi
  rm -rf "${EXPORTS:?}/$RUN"
fi
rm -f "$EXPORTS.$RUN.log"
# and the node trace (the tree viewer's input) with its records when it is small enough for git,
# else it stays in the store and the export. From the export when it has it (the same lines from the
# same units: a second read of every unit took tai60b's collect 74 minutes, then timed out); else
# qapctl nodes, which needs no --trace export
TRACED=0; [ -f "$EXPORTS/$RUN/nodes.txt.gz" ] && TRACED=1
# the node trace once the operation is over: a finished solve (results.txt), or one paused for an
# estimate or stopped (results_partial.txt)
if [ -f "$L/results.txt" ] || [ -f "$L/results_partial.txt" ]; then
  if [ "$TRACED" = 1 ]; then
    if [ "$(wc -c < "$EXPORTS/$RUN/nodes.txt.gz")" -gt 25000000 ]; then
      echo "node trace: over 25 MB gzipped, left in the export ($EXPORTS/$RUN/nodes.txt.gz) and the store"
    else
      cp "$EXPORTS/$RUN/nodes.txt.gz" "$L/nodes.txt.gz" && echo "node trace: from the export, $(du -h "$L/nodes.txt.gz" | cut -f1)"
    fi
  else
    NLOG=$(mktemp)
    if (cd "$D/../.." && python3 -m coord.qapctl nodes "$RUN" --backend "$BKS" \
          --out "$L/nodes.txt.gz") > "$NLOG" 2>&1; then
      if [ "$(wc -c < "$L/nodes.txt.gz")" -gt 25000000 ]; then
        rm -f "$L/nodes.txt.gz"
        echo "node trace: over 25 MB gzipped, left in the store (qapctl nodes $RUN --out FILE)"
      else
        echo "node trace: $(tail -1 "$NLOG" | sed 's/^[0-9:]* //')"
      fi
    else
      rm -f "$L/nodes.txt.gz" "$L/nodes.txt.gz.tmp"
      if grep -q "no frontier/nodes.txt" "$NLOG"; then :
      else echo "node trace: $(tail -1 "$NLOG")"; MISSING+=("node trace"); fi
    fi
    rm -f "$NLOG"
  fi
fi
get compare.txt && { echo; cat "$L/compare.txt"; }
echo
echo "logs: $L/hosts/*/out/ (node.log, coordinator.log, worker.log, gate_*.log)"
# the readout, made now (part-way included: an estimate's pending units imputed), into the record
echo
if (cd "$D/../.." && python3 -m coord.qapctl readout "$RUN" --backend "$BKS" --save \
      --out "$L/readout.md") > /dev/null 2>&1 && [ -s "$L/readout.md" ]; then
  cat "$L/readout.md"
else
  MISSING+=("readout"); echo "readout: could not be made (python3 -m coord.qapctl readout $RUN ...)"
fi
VERDICT=$(sed -n '3p' "$L/readout.md" 2> /dev/null | sed 's/\*\*//g; s/:.*//')
[ -n "$VERDICT" ] || VERDICT=UNKNOWN
if [ ${#MISSING[@]} -eq 0 ]; then C="collect complete"; else C="collect INCOMPLETE: no ${MISSING[*]}"; fi
echo
echo "$(date -u +%FT%TZ) $RUN: $VERDICT; $C" | tee "$L/collect.txt"
