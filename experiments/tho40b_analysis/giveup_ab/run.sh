#!/bin/bash
# The DNN give-up A/B on one machine (README.md here): the same complete subtree,
# searched three times with tho40-b's settings (infra/phase1/recipes/tho40.par) but
# for the give-up, through the coordinator's local mode (qapctl local), one arm
# after another:
#   g200    dnn.giveup = 200, dnn.gu_zone = 0.003   tho40-b's settings
#   g100    dnn.giveup = 100, dnn.gu_zone = 0.003   give up at 100; within the zone, at 200
#   g100z0  dnn.giveup = 100, dnn.gu_zone = 0       give up at 100, no zone
# (any arm gN or gNz0: dnn.giveup = N, with the zone or without it)
#   cut06   the grid cut phase at zone 0.06, the bound_lab dives' arm (no patch)
#   np10    bound_lab's ten-line block, with dnn.gu_zone = 0 and wc = 2.88 (no patch)
#   full12  that block with cut.handover and cut.cp_wc_gap: the integration build b86352e
# then prints a table of the arms' nodes, fathoms and CPU against g200's.
#
# From the repository's root, on the laptop (about 1 h an arm on 12 cores):
#   experiments/tho40b_analysis/giveup_ab/run.sh
# Options (environment): SUB (the -F subtree; default 26:10,8:0, the tho40ab subtree:
# 22,650 nodes, 78,804 c7a CPU s with -u 240516), UB (240516), SLOTS (all cores), ARMS
# ("g200 g100 g100z0"), OUT (runs/ here). An arm whose results are there is skipped,
# so the script can be run again after an interruption; a half-finished arm starts over.
set -euo pipefail
cd "$(dirname "$0")/../../.."
R=$(pwd)
SUB=${SUB:-26:10,8:0}
UB=${UB:-240516}
SLOTS=${SLOTS:-$( (sysctl -n hw.ncpu 2>/dev/null || nproc) )}
ARMS=${ARMS:-"g200 g100 g100z0"}
OUT=${OUT:-$R/experiments/tho40b_analysis/giveup_ab/runs}
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python3}
mkdir -p "$OUT/params"
say() { echo "$(date +%H:%M:%S) $*" | tee -a "$OUT/run.log"; }

say "give-up A/B: subtree -F $SUB, -u $UB, $SLOTS slots, arms: $ARMS (commit $(git rev-parse --short HEAD))"
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { say "build failed: $OUT/build.log"; exit 1; }

# the arms' params: tho40-b's, with the give-up changed
"$PY" - "$R/infra/phase1/recipes/tho40.par" "$OUT/params" $ARMS << 'EOF'
import re, sys, os
base, out, names = sys.argv[1], sys.argv[2], sys.argv[3:]
# the grid cut phase arms of experiments/bound_lab (README.md, section 2; tho40_dives/par/):
# cut06 the dives' zone 0.06 arm (226 c7a CPU-days by dives), no patch; np10 the README's
# ten-line block that needs no patch, with dnn.gu_zone = 0 and the phase penalty; full12 that
# block with the two keys of the integration build (b86352e), the dives' combined arm
CUT06 = {"pool.on": "1", "cut.grid": "5:8", "cut.side": "1", "cut.cp": "100:1000", "cut.warmy": "1",
         "cut.rounds": "12", "cut.iters": "200", "cut.freeze": "1000", "cut.zone": "0.06"}
NP10 = dict(CUT06, **{"dnn.giveup_warm": "100", "dnn.gu_zone": "0", "dnn.beta_sched": "w=0.36,wc=2.88"})
FULL12 = dict(NP10, **{"cut.handover": "20:0:0.22:1.5:100:3:1.25", "cut.cp_wc_gap": "300"})
NAMED = {"cut06": CUT06, "np10": NP10, "full12": FULL12}
arms = {}
for a in names:   # gN: dnn.giveup = N with the zone (0.003); gNz0: without it; or a named arm
    m = re.fullmatch(r"g(\d+)(z0)?", a)
    if a in NAMED:
        arms[a] = NAMED[a]
    elif m:
        arms[a] = {"dnn.giveup": m.group(1), "dnn.gu_zone": "0" if m.group(2) else "0.003"}
    else:
        sys.exit(f"arm {a}: not gN, gNz0, cut06, np10 or full12")
text = open(base).read()
for arm, kv in arms.items():
    t = text
    for k, v in kv.items():
        pat = rf"^{re.escape(k)} = .*$"
        t = re.sub(pat, f"{k} = {v}", t, flags=re.M) if re.search(pat, t, re.M) else t.rstrip("\n") + f"\n{k} = {v}\n"
    open(os.path.join(out, arm + ".par"), "w").write(t)
EOF

# full12 needs the integration build (experiments/bound_lab/integration/README.md): its own
# worktree at INTEGRATION, built once; the coordinator (coord/) is this checkout's, unchanged there
INTEGRATION=${INTEGRATION:-b86352e}
INTSOLVER=""
case " $ARMS " in *" full12 "*)
  WT=${INTEGRATION_WT:-$R/../qap-$INTEGRATION}
  if [ ! -x "$WT/latest/bnbdrv" ]; then
    git fetch -q origin claude/cool-babbage-ilnogj-integration || true
    [ -d "$WT" ] || git worktree add --detach "$WT" "$INTEGRATION" >> "$OUT/build.log" 2>&1
    make -C "$WT" solver PYTHON="$PY" >> "$OUT/build.log" 2>&1 || { say "integration build failed: $OUT/build.log"; exit 1; }
  fi
  INTSOLVER="$WT/latest/bnbdrv"
  say "full12: solver $INTSOLVER ($(git -C "$WT" rev-parse --short HEAD))" ;;
esac

for arm in $ARMS; do
  res="$OUT/$arm.results.txt"
  solver="$R/latest/bnbdrv"; [ "$arm" = full12 ] && solver="$INTSOLVER"
  if [ -s "$res" ]; then say "$arm: done before, skipped"; continue; fi
  d="$OUT/$arm"; rm -rf "$d"
  say "$arm: running ($(grep -E '^(dnn\.(giveup|gu_zone|giveup_warm|beta_sched)|cut\.(zone|handover|cp_wc_gap)) ' "$OUT/params/$arm.par" | tr '\n' ' '))"
  if "$PY" -m coord.qapctl local "gu-$arm" --backend "local:$d" --instance "$R/qaplib/tho40t.dat" \
       --params "$OUT/params/$arm.par" --ub "$UB" --slots "$SLOTS" --solver "$solver" \
       --extra "-F $SUB -s $R/latest/SYM/tho40t.sym" \
       --budget "${BUDGET:-300}" --target-unit-cpu "${BUDGET:-300}" --steal-after 60 --checkpoint-every 60 --ckpt-upload 120 \
       --trace --visibility private --workdir "$d/work" > "$OUT/$arm.log" 2>&1; then
    cp "$d/store/runs/gu-$arm/results.txt" "$res"
    say "$arm: $(grep -E 'Nodes:|coordinator:' "$res" | tr -s ' ' | tr '\n' ' ' | cut -c1-200)"
    # the node and DNN traces (coord/export.py), then the run's store (attempts, states): its results are copied
    "$PY" -m coord.qapctl export "gu-$arm" --backend "local:$d" --out "$OUT/$arm.export" >> "$OUT/$arm.log" 2>&1 \
      && rm -rf "$d" || say "$arm: export failed, the store is kept in $d"
  else
    say "$arm: FAILED, see $OUT/$arm.log"
  fi
done

"$PY" - "$OUT" $ARMS << 'EOF' | tee "$OUT/summary.txt"
import re, sys, os
out, arms = sys.argv[1], sys.argv[2:]
rows = {}
for a in arms:
    p = os.path.join(out, a + ".results.txt")
    if not os.path.exists(p):
        continue
    t = open(p).read()
    n = re.search(r"Nodes:\s*(\d+)\s*/\s*(\d+)", t)
    c = re.search(r"wall ([\d.]+) s, CPU ([\d.]+) s", t)
    rows[a] = (int(n.group(1)), int(n.group(2)), float(c.group(1)), float(c.group(2)))
print(f"{'arm':>7} {'nodes':>8} {'fathomed':>9} {'wall s':>8} {'CPU s':>10} {'nodes / g200':>13} {'CPU / g200':>11}")
b = rows.get("g200")
for a, (n, f, w, c) in rows.items():
    print(f"{a:>7} {n:>8,} {f:>9,} {w:>8,.0f} {c:>10,.0f} {n / b[0] if b else float('nan'):>13.3f} {c / b[3] if b else float('nan'):>11.3f}")
EOF
say "summary: $OUT/summary.txt"
