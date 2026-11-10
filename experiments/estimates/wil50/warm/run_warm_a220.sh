#!/bin/bash
# a220's whole subtree, warm: its one-node state cut from r.12.22's out.full.state (the
# parent bounded as production does under fullrun.par, its phase and pool included), then
# parallel.py --resume on 4 cores with fullrun.par. Resumable: run again (checkpoint).
S=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/warm220
R=/home/user/qap
cd $R || exit 1
. experiments/wil50_feasibility/fullrun/settings.sh || exit 1
fullrun_check || exit 1
if [ ! -f $S/node.state ]; then
  python3 - $S/stageA/tasks/r.12.22/out.full.state $S/node.state <<'PY' || exit 1
import sys
sys.path.insert(0, "experiments/wil50_feasibility/sampling")
import common
st = common.State(sys.argv[1])
idx = [i for i, (t, w) in enumerate(st.nodes) if common.spec(common.node_fields(t, st.n)[3]) == "21:22,31:23,40:15"]
if len(idx) != 1:
    raise SystemExit(f"a220 is not among r.12.22's children here ({len(st.nodes)} children)")
print("a220 is child", idx[0], "of", len(st.nodes), common.node_fields(st.nodes[idx[0]][0], st.n)[:3])
common.write_one(sys.argv[1], idx[0], sys.argv[2])
PY
fi
# the state records stage A copy of fullrun.par (same settings, a different par.hash)
RES=(); [ -f $S/warm.ckpt ] && RES=(--resume $S/warm.ckpt) || RES=(--resume $S/node.state)
QAP_DNN_TRACE=$S/trace.txt python3 bench/parallel.py "$FULLRUN_INSTANCE" --jobs 4 --ub "$FULLRUN_UB" --params "$FULLRUN_PAR" \
  --extra "-s $FULLRUN_SYM --allow-settings-change" "${RES[@]}" --checkpoint $S/warm.ckpt --checkpoint-every 900 --wall-limit ${WALL:-6600} \
  --workdir $S/work >> $S/parallel.log 2>&1
echo "parallel exit $?" | tee -a $S/parallel.log
tail -3 $S/parallel.log
