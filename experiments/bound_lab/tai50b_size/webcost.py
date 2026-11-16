"""E3, the web model's tree priced with measured costs [e]:  python3 webcost.py  -> logs/webcost.txt
nodes(d), open(d) from logs/web.txt (web.py) for the gap scales 0.5 and 1.0, depths 1..8;
an open node costs what the spine node of its depth cost [m]; a fathomed node the mean CPU of the fathomed node
bounds of that depth in the dive sets [m] (depth 1: of the sampled depth-1 nodes that fathomed);
below each open depth-8 node: LA_tail8's mean dive estimate; plus the root, and the children of the depth-1
near misses the model does not have (4 of the 15 sampled non-twin nodes branched at the phase limit or the stall
rule: 45 x 4/15 = 12 parents x 49 children at the mean cost measured under them, sets LB2 and LB_near)."""
import os, re, sys, glob, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, dives as DV, estimate as E
H = os.path.dirname(os.path.abspath(__file__)); cpu = lambda d: d["fcpu"] + d["ccpu"]
root = N.load([H + "/runs/root_long/dnn.txt"])[0]
d1 = [d for d in N.load(sorted(glob.glob(H + "/runs/L1_*/dnn.txt"))) if d["depth"] == 1]
sp = {d["depth"]: d for d in N.load([H + "/runs/L1_0/dnn.txt"] + sorted(glob.glob(H + "/runs/spine_*/dnn.txt"))) if d["onpath"]}
allb = N.load([f for f in glob.glob(H + "/runs/L*/s*/dnn.txt")])
cf = {1: st.mean(cpu(d) for d in d1 if d["fath"])}; nf = {1: sum(d["fath"] for d in d1)}
for dep in range(2, 9):
    X = [cpu(d) for d in allb if d["depth"] == dep and d["fath"]]; cf[dep] = st.mean(X) if X else float("nan"); nf[dep] = len(X)
tail = E.dive_rows("LA_tail8"); tn = st.mean(r[1] for r in tail) if tail else float("nan"); tc = st.mean(r[2] for r in tail) if tail else float("nan")
near = [cpu(d) for s in ("LB2", "LB_near") for d in N.load(glob.glob(f"{H}/runs/{s}/s*/dnn.txt")) if d["depth"] == 2]
out = ["cost of an open node by depth (the spine's) [m]: " + ", ".join("%d: %.0f" % (k, cpu(v)) for k, v in sorted(sp.items())),
       "mean cost of a fathomed node by depth [m] (node bounds): " + ", ".join("%d: %.0f (%d)" % (k, cf[k], nf[k]) for k in sorted(cf)),
       "below an open depth-8 node (LA_tail8, %d dives) [e]: %.0f nodes, %.0f CPU-s" % (len(tail), tn, tc),
       "children of the depth-1 near misses: 12 parents x 49 x %.0f CPU-s (%d node bounds) [e]" % (st.mean(near), len(near))]
blocks = re.split(r"\ngap scale ", open(H + "/logs/web.txt").read())[1:]
for b in blocks:
    scale = float(b.split(":")[0]); rows = [[int(x) for x in ln.split("|")[1:4]] for ln in b.splitlines() if re.match(r"\| \d", ln)]
    if "stopped" in b: out.append("\ngap scale %.1f: the model was stopped at depth %d with %d open nodes" % (scale, rows[-1][0], rows[-1][2])); continue
    nodes = 1 + 12 * 49; c = cpu(root) + 12 * 49 * st.mean(near); tab = []
    for dep, n, op in rows:
        if dep > 8: continue
        cc = op * cpu(sp[dep]) + (n - op) * cf[dep]; nodes += n; c += cc; tab.append("| %d | %d | %d | %.0f |" % (dep, n, op, cc))
    o8 = [r for r in rows if r[0] == 8][0][2]; nodes += o8 * tn; c += o8 * tc
    out += ["\ngap scale %.1f: **%.0f nodes, %.0f CPU-s = %.0f core-hours**" % (scale, nodes, c, c / 3600), "", "| depth | nodes | open | CPU-s |", "|---|---|---|---|",
            "| 0 | 1 | 1 | %.0f |" % cpu(root)] + tab + ["| 2, under the 12 near misses | %d | 0 | %.0f |" % (12 * 49, 12 * 49 * st.mean(near)), "| 9 and below | %.0f | | %.0f |" % (o8 * tn, o8 * tc)]
open(H + "/logs/webcost.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
