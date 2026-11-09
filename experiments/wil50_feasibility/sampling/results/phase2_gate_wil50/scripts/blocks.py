"""which stage-A nodes would change parameter block (pick_param) if the root bound were 18ebd58's"""
import json, os
SA = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA/tasks"
opt = 48816.0
R_OLD, R_NEW = 46100.957091528995, 46131.507862715953
# blocks: (relative_gap, depth) of fullrun.par / nug_deep_xbar9_all.par; the last is the bottom rule
BLOCKS = [(0.20, 2), (0.12, 4), (0.0, 6), (-1, 1000)]


def block(bnd, depth, root):
    g = (opt - bnd) / (opt - root)
    i = 0
    while i < len(BLOCKS) - 1 and g <= BLOCKS[i][0]:
        i += 1
    while i < len(BLOCKS) - 1 and depth > BLOCKS[i][1]:
        i += 1
    return i


moved = {}
for t in os.listdir(SA):
    d = json.load(open(os.path.join(SA, t, "done.json")))
    kids = d.get("children") or d.get("frontier") or []
    for k in kids:
        depth = d["depth"] + 1
        b = max(k["bound"], k["pred"])
        a, c = block(b, depth, R_OLD), block(b, depth, R_NEW)
        if a != c:
            moved.setdefault((d["depth"] + 1, a, c), []).append((t, round(k["bound"], 4)))
for key, v in sorted(moved.items()):
    parents = sorted({p for p, _ in v})
    print(f"depth {key[0]}: block {key[1]} -> {key[2]}: {len(v)} nodes, parents {parents[:12]}{' ...' if len(parents) > 12 else ''} bounds {sorted({b for _, b in v})[:6]}")
print("thresholds old/new: 0.2 ->", opt - 0.2 * (opt - R_OLD), opt - 0.2 * (opt - R_NEW), "; 0.12 ->",
      opt - 0.12 * (opt - R_OLD), opt - 0.12 * (opt - R_NEW))
