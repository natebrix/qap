import json, os, sys, collections
S = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA/tasks"
rows = []
for t in os.listdir(S):
    p = os.path.join(S, t, "done.json")
    d = json.load(open(p))
    if d["depth"] != 2:
        continue
    rows.append((t, d, d["trace"]))
print(len(rows), "depth-2 tasks")
c = collections.Counter()
for t, d, tr in rows:
    if not tr:
        key = ("noDNN",)
    else:
        x = tr[0]
        key = (x["status"], x["cut"], x.get("cut_status"), x["fathomed"], len(tr))
    c[key] += 1
for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
    print(v, k)
