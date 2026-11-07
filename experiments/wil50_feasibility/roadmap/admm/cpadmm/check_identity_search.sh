#!/bin/bash
# The search with the CP switches off is unchanged: this tree's bnbdrv against the
# parent commit's (7849e38, a scratch build; OLD_BNB), per-node trace lines
# (QAP_DNN_TRACE, pid and CPU stripped) and the search's node counts:
#   nug20   the ship test's DNN settings (cut phases at zone 0.005 on the 4 x 5 grid,
#           pools, the pool check), the whole search
#   node10  the bestcore arm's environment (pipeline_ab/run_pair.sh) on the node-10
#           subtree (-F 39:27,30:48,25:49), a CPU slice of SLICE seconds (default 900):
#           the nodes both traces have must be identical, and both must have the
#           same first nodes (the slice ends at a different node only by timing)
# Cores: CORE_OLD, CORE_NEW (0, 1). Logs in logs/ident_search/.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
OLD_BNB=${OLD_BNB:-/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/bnbdrv_old}
NEW_BNB=$R/latest/bnbdrv
SLICE=${SLICE:-900}
C0=${CORE_OLD:-0}; C1=${CORE_NEW:-1}
L=$D/logs/ident_search; mkdir -p "$L"
cd "$R/latest" || exit 1
SHIP_ENV="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 DNN_CHOLQR_PASSES=1 QAP_DNN_WARM_Y=2 QAP_DNN_EXT=20 QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=4:5 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20"
# run_pair.sh's arm_env base + bestcore
BEST_ENV="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20 QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1000 DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 QAP_DNN_WARM_Y=2 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200 QAP_DNN_CUT_DAMP1=0.3 QAP_DNN_SHIP=32 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
strip() { grep ' fixes ' "$1" | sed -E 's/^pid [0-9]+ //; s/ cpu [0-9.]+//g'; }
# nug20, the whole search
for v in old new; do
  bin=$OLD_BNB; c=$C0; [ $v = new ] && { bin=$NEW_BNB; c=$C1; }
  rm -f "$L/nug20_$v.trace"
  env $SHIP_ENV QAP_DNN_TRACE="$L/nug20_$v.trace" taskset -c $c "$bin" -p PARAM/nug_deep_xbar9_all.par -x 0:4:1:300:1 -u 2570 \
    ../qaplib/nug20.dat > "$L/nug20_$v.out" 2>&1 &
done
wait
if cmp -s <(strip "$L/nug20_old.trace") <(strip "$L/nug20_new.trace") && \
   [ "$(grep 'Nodes:' "$L/nug20_old.out")" = "$(grep 'Nodes:' "$L/nug20_new.out")" ]; then
  echo "nug20 (ship settings): identical ($(strip "$L/nug20_old.trace" | wc -l) node lines; $(grep 'Nodes:' "$L/nug20_new.out" | sed 's/^ *//'))"
else echo "nug20: DIFFERS"; diff <(strip "$L/nug20_old.trace") <(strip "$L/nug20_new.trace") | head -6; fi
# node 10, a CPU slice of the bestcore arm
for v in old new; do
  bin=$OLD_BNB; c=$C0; [ $v = new ] && { bin=$NEW_BNB; c=$C1; }
  rm -f "$L/node10_$v.trace" "$L/node10_$v.state"
  env $BEST_ENV QAP_DNN_TRACE="$L/node10_$v.trace" taskset -c $c "$bin" -p PARAM/nug_deep_xbar9_all.par -u 48816 \
    -F 39:27,30:48,25:49 -x 0:4:1:300:1 -o "$L/node10_$v.state" -t $SLICE ../qaplib/wil50.dat > "$L/node10_$v.out" 2>&1 &
done
wait
rm -f "$L"/node10_*.state
python3 - "$L/node10_old.trace" "$L/node10_new.trace" <<'EOF'
import re, sys
def rows(p):
    out = []
    for l in open(p):
        if " fixes " in l:
            out.append(re.sub(r" cpu \S+", "", re.sub(r"^pid \d+ ", "", l)).strip())
    return out
a, b = rows(sys.argv[1]), rows(sys.argv[2])
k = min(len(a), len(b))
same = a[:k] == b[:k]
print("node10 (bestcore settings, a %s-node and a %s-node slice): the first %d node lines %s; phases among them %d"
      % (len(a), len(b), k, "identical" if same else "DIFFER", sum(" | cut 1" in x for x in a[:k])))
if not same:
    for x, y in zip(a, b):
        if x != y: print("  old:", x[:200]); print("  new:", y[:200]); break
EOF
