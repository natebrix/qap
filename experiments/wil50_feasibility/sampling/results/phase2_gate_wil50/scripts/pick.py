import json, os
S = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA/tasks"
cls = {"first_fathomed": [], "cp_fathomed": [], "stalled": [], "iter_limit": [], "hopeless_nocut": []}
for t in os.listdir(S):
    d = json.load(open(os.path.join(S, t, "done.json")))
    if d["depth"] != 2:
        continue
    x = d["trace"][0]
    if x["status"] == "fathomed" and not x["cut"]:
        k = "first_fathomed"
    elif x["cut"] and x.get("cut_status") == "fathomed":
        k = "cp_fathomed"
    elif x["cut"] and x.get("cut_status") == "stalled":
        k = "stalled"
    elif x["cut"] and x.get("cut_status") == "iter_limit":
        k = "iter_limit"
    else:
        k = "hopeless_nocut"
    cls[k].append((d["stats"]["times"][0], t, d["spec"], x["iters"], x.get("cut_iters"), x.get("cut_cp"),
                   round(x["lb"], 4), round(x.get("cut_best", 0), 4), d["n_children"], d["line"], d.get("retained")))
for k, v in cls.items():
    v.sort()
    print(k, len(v))
    for r in v[:6]:
        print("   cpu %.1f  %s  [%s]  it %s cut_it %s cp %s lb %s best %s kids %s %s retained %s" % r)
