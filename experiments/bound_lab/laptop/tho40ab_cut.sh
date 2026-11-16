#!/bin/bash
# The cut-phase arms on the tho40ab subtree, through the coordinator's local mode, to sit
# next to experiments/tho40b_analysis/giveup_ab (same subtree, same command, same table):
#   cut06   tho40-b's settings + the grid cut phase at cut.zone = 0.06 (settings only)
#   comb    cut06 + early entry, the hand-over rule and the phase's own penalty
#           (bound_lab/README.md section 2; needs this branch's solver: cut.handover, cut.cp_wc_gap)
#   g200    tho40-b's settings (the reference; skipped when giveup_ab's result is there)
# From the repository's root, on the laptop:
#   experiments/bound_lab/laptop/tho40ab_cut.sh
# Options (environment): SUB (default 26:10,8:0: 22,650 nodes, 78,804 c7a CPU s under tho40-b's
# settings), UB (240516), SLOTS (all cores), ARMS ("cut06 comb"), OUT (runs_tho40ab/ here),
# REF (a results.txt of the g200 arm; default giveup_ab's in this checkout or in $COGNAC).
# An arm whose results are there is skipped; a half-finished arm starts over.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../.."
R=$(pwd)
SUB=${SUB:-26:10,8:0}
UB=${UB:-240516}
SLOTS=${SLOTS:-$( (sysctl -n hw.ncpu 2>/dev/null || nproc) )}
ARMS=${ARMS:-"cut06 comb"}
OUT=${OUT:-$H/runs_tho40ab}
REF=${REF:-}
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python3}
mkdir -p "$OUT/params"
say() { echo "$(date +%H:%M:%S) $*" | tee -a "$OUT/run.log"; }

say "tho40ab cut arms: subtree -F $SUB, -u $UB, $SLOTS slots, arms: $ARMS (commit $(git rev-parse --short HEAD))"
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { say "build failed: $OUT/build.log"; exit 1; }

"$PY" - "$R/infra/phase1/recipes/tho40.par" "$OUT/params" << 'EOF'
import re, sys, os
base, out = sys.argv[1], sys.argv[2]
cut = {"pool.on": "1", "cut.grid": "5:8", "cut.side": "1", "cut.cp": "100:1000", "cut.warmy": "1",
       "cut.rounds": "12", "cut.iters": "200", "cut.freeze": "1000", "cut.zone": "0.06"}
comb = dict(cut, **{"dnn.gu_zone": "0", "dnn.giveup_warm": "100", "dnn.beta_sched": "w=0.36,wc=2.88",
                    "cut.handover": "20:0:0.22:1.5:100:3:1.25", "cut.cp_wc_gap": "300"})
arms = {"g200": {}, "cut06": cut, "comb": comb}
text = open(base).read()
for arm, kv in arms.items():
    t = text
    for k, v in kv.items():
        pat = rf"^{re.escape(k)} = .*$"
        t = re.sub(pat, f"{k} = {v}", t, flags=re.M) if re.search(pat, t, re.M) else t.rstrip("\n") + f"\n{k} = {v}\n"
    open(os.path.join(out, arm + ".par"), "w").write(t)
EOF

for arm in $ARMS; do
  res="$OUT/$arm.results.txt"
  if [ -s "$res" ]; then say "$arm: done before, skipped"; continue; fi
  d="$OUT/$arm"; rm -rf "$d"
  say "$arm: running"
  if "$PY" -m coord.qapctl local "cut-$arm" --backend "local:$d" --instance "$R/qaplib/tho40t.dat" \
       --params "$OUT/params/$arm.par" --ub "$UB" --slots "$SLOTS" \
       --extra "-F $SUB -s $R/latest/SYM/tho40t.sym" \
       --budget 300 --target-unit-cpu 300 --steal-after 60 --checkpoint-every 60 --ckpt-upload 120 \
       --trace --visibility private --workdir "$d/work" > "$OUT/$arm.log" 2>&1; then
    cp "$d/store/runs/cut-$arm/results.txt" "$res"
    say "$arm: $(grep -E 'Nodes:|coordinator:' "$res" | tr -s ' ' | tr '\n' ' ' | cut -c1-200)"
    "$PY" -m coord.qapctl export "cut-$arm" --backend "local:$d" --out "$OUT/$arm.export" >> "$OUT/$arm.log" 2>&1 \
      && rm -rf "$d" || say "$arm: export failed, the store is kept in $d"
  else
    say "$arm: FAILED, see $OUT/$arm.log"
  fi
done

# the reference: this script's g200 if it ran, else giveup_ab's
if [ -z "$REF" ]; then
  for c in "$OUT/g200.results.txt" "$R/experiments/tho40b_analysis/giveup_ab/runs/g200.results.txt" \
           "${COGNAC:-$R/../qap}/experiments/tho40b_analysis/giveup_ab/runs/g200.results.txt"; do
    [ -s "$c" ] && { REF=$c; break; }
  done
fi

"$PY" - "$OUT" "${REF:-none}" g200 $ARMS << 'EOF' | tee "$OUT/summary.txt"
import re, sys, os
out, ref, arms = sys.argv[1], sys.argv[2], list(dict.fromkeys(sys.argv[3:]))
def read(p):
    t = open(p).read()
    n = re.search(r"Nodes:\s*(\d+)\s*/\s*(\d+)", t)
    c = re.search(r"wall ([\d.]+) s, CPU ([\d.]+) s", t)
    return (int(n.group(1)), int(n.group(2)), float(c.group(1)), float(c.group(2)))
rows = {}
for a in arms:
    p = os.path.join(out, a + ".results.txt")
    if a == "g200" and not os.path.exists(p) and ref != "none":
        p = ref
    if os.path.exists(p):
        rows[a] = read(p)
b = rows.get("g200")
print(f"reference (g200): {ref if b else 'none found: run with ARMS=\"g200 cut06 comb\" or set REF'}")
print(f"{'arm':>7} {'nodes':>8} {'fathomed':>9} {'wall s':>8} {'CPU s':>10} {'g200 nodes / arm':>17} {'g200 CPU / arm':>15}")
for a, (n, f, w, c) in rows.items():
    print(f"{a:>7} {n:>8,} {f:>9,} {w:>8,.0f} {c:>10,.0f} {b[0] / n if b else float('nan'):>17.2f} {b[3] / c if b else float('nan'):>15.2f}")
EOF
say "summary: $OUT/summary.txt"
