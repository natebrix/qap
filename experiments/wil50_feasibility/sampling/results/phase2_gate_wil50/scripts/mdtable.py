import json, os, sys
W = "/home/user/qap/.claude/worktrees/agent-aebd90129de6196e4"
SA = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA"
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
sys.path.insert(0, os.path.join(W, "experiments", "wil50_feasibility", "sampling"))
from common import State, node_fields  # noqa: E402

rows = {(r["side"], r["task"]): r for r in json.load(open(G + "/compare.json"))}
order = ["r.9.39", "r.4.39", "r.3.36", "r.7.35", "r.6.41", "r.2.26", "r.3.15", "r.13.4", "r.1", "r"]
CLS = {"r.9.39": "first solve", "r.4.39": "first solve", "r.3.36": "CP fathomed", "r.7.35": "CP fathomed",
       "r.6.41": "stalled", "r.2.26": "stalled", "r.3.15": "iter limit", "r.13.4": "iter limit", "r.1": "depth 1",
       "r": "root, QPB"}


def old_pred(tid):
    for name in ("out.full.state", "out.state", "out.slim.state"):
        p = os.path.join(SA, "tasks", tid, name)
        if os.path.exists(p):
            st = State(p)
            pr = [node_fields(t, st.n)[2] for t, _ in st.nodes]
            return (min(pr), max(pr), st.root_bound) if pr else None
    return None


def fmt(t, r, kids_bound, pred, cpu, line, nk):
    if t:
        first = f"{t['lb']:.4f}"
        st = f"{t['status']} {t['iters']}"
        cut = f"{t['cut_status']} {t['cut_best']:.4f}, {t['cut_iters']} it, {t['cut_rounds']} rd" if t["cut"] else "-"
        fat = str(t["fathomed"])
    else:
        first, st, cut, fat = "-", "QPB only", "-", "0"
    kb = f"{kids_bound!r}" if kids_bound is not None else "-"
    pr = f"{pred[0]:.2f} - {pred[1]:.2f}" if pred else "-"
    return f"{first} | {st} | {cut} | {fat} | {nk}, {line} | {kb} | {pr} | {cpu:.1f}"


print("| task | class | side | first lb | status, iters | cut phase: outcome, best, iters, rounds | fathomed | children, line | children's inherited bound | children's predicted bounds | CPU s |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for tid in order:
    r = rows[("new", tid)]
    op = old_pred(tid)
    ob = None
    if r["n_children_record"]:
        # the record's inherited bound: from its own file
        for name in ("out.full.state", "out.state", "out.slim.state"):
            p = os.path.join(SA, "tasks", tid, name)
            if os.path.exists(p):
                st = State(p)
                ob = sorted({node_fields(t, st.n)[1] for t, _ in st.nodes})
                break
    old = fmt(r["rec_trace"], r, (ob[0] if ob and len(ob) == 1 else ob), op[:2] if op else None, r["cpu_record"],
              r["line_record"], r["n_children_record"])
    nb = r.get("children_bound")
    new = fmt(r["rep_trace"], r, (nb[0] if nb and len(nb) == 1 else nb), r.get("children_pred_min_max"), r["cpu_replay"],
              r["line_replay"], r["n_children_replay"])
    print(f"| `{tid}` `{r['spec'] or 'root'}` | {CLS[tid]} | stage A | {old} |")
    print(f"| | | HEAD | {new} |")
    if tid == "r":
        print(f"<!-- root bound record {op[2]} -->")
