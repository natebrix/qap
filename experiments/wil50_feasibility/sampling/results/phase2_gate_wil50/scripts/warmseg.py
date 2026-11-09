import os, sys
W = "/home/user/qap/.claude/worktrees/agent-aebd90129de6196e4"
SA = "/home/user/qap/.claude/worktrees/agent-abe2ac37e80fb6b2b/experiments/wil50_feasibility/sampling/runs/stageA"
G = "/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/gate"
sys.path.insert(0, os.path.join(W, "experiments", "wil50_feasibility", "sampling"))
from common import State  # noqa: E402

for side, tid in [a.split(":") for a in sys.argv[1:]]:
    old_p = os.path.join(SA, "tasks", tid, "out.full.state")
    if not os.path.exists(old_p):
        old_p = os.path.join(SA, "tasks", tid, "out.state")
    so, sn = State(old_p), State(os.path.join(G, "runs", side, tid, "out.state"))
    wo = {id(w): w for _, w in so.nodes if w}
    wn = {id(w): w for _, w in sn.nodes if w}
    assert len(wo) == 1 and len(wn) == 1, (len(wo), len(wn))
    a, b = list(wo.values())[0].text.split(" "), list(wn.values())[0].text.split(" ")
    m = int(a[0]); n = so.n
    segs = {"m": (0, 1), "p": (1, 1 + n), "X": (1 + n, 1 + n + m * m), "dnn_ext": (1 + n + m * m, 2 + n + m * m),
            "pool": (2 + n + m * m, 3 + n + m * m), "dnn_state": (3 + n + m * m, 4 + n + m * m)}
    print(f"{side} {tid}: {len(a)}/{len(b)} fields, m {m}")
    for k, (i, j) in segs.items():
        x, y = a[i:j], b[i:j]
        same = x == y
        extra = ""
        if k == "X" and not same:
            fx, fy = [float(v) for v in x], [float(v) for v in y]
            extra = f" max |dX| {max(abs(p - q) for p, q in zip(fx, fy)):.3g}, {sum(p != q for p, q in zip(fx, fy))} of {len(fx)} entries differ"
        if k in ("pool", "dnn_state"):
            extra = f" ({len(x[0])} / {len(y[0])} chars)"
        print(f"   {k:9s} {'same' if same else 'DIFF'}{extra}")
    # the node lines: which fields differ
    nb = sum(1 for (t1, _), (t2, _) in zip(so.nodes, sn.nodes) if t1.split()[2] != t2.split()[2])
    npd = sum(1 for (t1, _), (t2, _) in zip(so.nodes, sn.nodes) if t1.split()[3] != t2.split()[3])
    nrest = sum(1 for (t1, _), (t2, _) in zip(so.nodes, sn.nodes) if t1.split()[:2] + t1.split()[4:] != t2.split()[:2] + t2.split()[4:])
    dp = [float(t2.split()[3]) - float(t1.split()[3]) for (t1, _), (t2, _) in zip(so.nodes, sn.nodes)]
    print(f"   nodes: inherited bound differs on {nb}, predicted bound on {npd} (new - old: min {min(dp):.3f} max {max(dp):.3f}), other fields on {nrest} of {len(so.nodes)}")
    print(f"   header root_bound old {so.root_bound} new {sn.root_bound}")
