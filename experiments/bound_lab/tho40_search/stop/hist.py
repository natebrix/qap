"""The recorded phase histories (QAP_DNN_TRACE_HIST: the running certified best every 10 iterations)
and what followed each phase in its run.

load_phases() -> list of dicts: run, depth, m, fixes (trace order key), pooled, short0 (shortfall of
the solve before the phase), cutoff, hist [(it, gap)], status, its, left (gap at the end, 0 if fathomed),
after (iterations spent under the node's children when the phase missed, None if not run / not complete).
"""
import os, re, glob, sys
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, H)
from summarize import load

def load_phases(runs=None):
    out = []
    for d in sorted(glob.glob(os.path.join(H, "runs", "*_*"))):
        name = os.path.basename(d)
        if runs and name not in runs: continue
        p = os.path.join(d, "dnn.txt")
        if not os.path.exists(p): continue
        nodes, dnn = load(d)
        lines = [l for l in open(p) if " | first lb " in l]
        below = {}
        for x, e in zip(nodes, dnn):
            f = x["fixes"].split(",")
            for k in range(1, len(f)):
                v = below.setdefault(",".join(f[:k]), [0, 0]); v[0] += 1; v[1] += e["iters"] + e["cit"]
        for x, e, ln in zip(nodes, dnn, lines):
            if not e["cut"] or "| hist" not in ln: continue
            h = [(int(a), e["cutoff"] - float(b)) for a, b in (t.split(":") for t in ln.split("| hist")[1].split())]
            fath = e["cstatus"] == "fathomed"
            out.append(dict(run=name, depth=x["depth"], m=40 - x["depth"], fixes=x["fixes"], pooled=e["kept"] > 0,
                            short0=e["cutoff"] - e["lb"], cutoff=e["cutoff"], hist=h, status=e["cstatus"], its=e["cit"],
                            left=0.0 if fath else e["cutoff"] - e["cbest"], first=e["iters"],
                            after=None if fath else below.get(x["fixes"], [0, None])[1], kids=below.get(x["fixes"], [0, 0])[0],
                            path=name.startswith("path")))
    return out

if __name__ == "__main__":
    for p in load_phases():
        g = dict(p["hist"])
        def at(t): return g.get(t, float("nan"))
        print("%-11s d%d %s short %7.1f -> %-9s %4d it left %7.1f | gap at 200 %6.0f 400 %6.0f 600 %6.0f 1000 %6.0f 1500 %6.0f | children after: %s nodes %s it" % (
            p["run"], p["depth"], "pool" if p["pooled"] else "    ", p["short0"], p["status"], p["its"], p["left"], at(200), at(400), at(600), at(1000), at(1500), p["kids"] if p["after"] is not None else "-", p["after"] if p["after"] is not None else "-"))
