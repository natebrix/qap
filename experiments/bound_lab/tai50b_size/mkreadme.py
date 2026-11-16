"""README.md from README.tmpl.md: the @NAME@ marks are filled from logs/ (run ./mkreport.sh first) and from the
hand-written answer.md and notdone.md.   python3 mkreadme.py"""
import os, glob, sys, statistics as st, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, dives as DV, estimate as E
H = os.path.dirname(os.path.abspath(__file__)); rd = lambda f: open(os.path.join(H, f)).read().rstrip("\n")
def nodetable(files):
    out = ["| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |", "|---|---|---|---|---|---|---|---|---|"]
    for d in N.load(files):
        end = d["cstatus"].replace("stalled stalled", "stalled").replace("iter_limit", "limit")
        out.append("| `%s` | %s | %d | %s | %s | %s | %s | %s | %.0f |" % (",".join(d["fixes"]), "yes" if d["onpath"] else "", d["fit"], f"{d['cutoff'] - d['lb']:,.0f}",
                   d["cit"] or "", f"{d['cutoff'] - d['best']:,.0f}" if d["cut"] else "", end, {"B": "branched", "C": "fathomed (phase)", "D": "fathomed (first solve)"}[d["out"]], d["fcpu"] + d["ccpu"]))
    return "\n".join(out)
def divetable():
    out = ["| set | start nodes (depth) | processes | dives | node bounds: fathomed by the first solve / by the phase / branched | deepest dive | mean CPU-s per node bound | estimate of the set's subtrees: nodes | CPU-s |", "|---|---|---|---|---|---|---|---|---|"]
    sets = [f"LA_off{d}" for d in range(2, 9)] + ["LA_tail8"] + sorted(os.path.basename(p) for p in glob.glob(H + "/runs/LB*"))
    for s in sets:
        rows = E.dive_rows(s)
        if not rows: out.append(f"| `{s}` | | | 0 | | | | | |"); continue
        ch = [d for f, c in DV.dives(f"{H}/runs/{s}") for d in c]
        st0 = [ln for ln in open(glob.glob(f"{H}/runs/{s}/s*/est.out")[0]) if ln.startswith("stack =")]
        se = lambda k: st.stdev(r[k] for r in rows) / math.sqrt(len(rows)) if len(rows) > 1 else float("nan")
        out.append("| `%s` | %s (%d) | %d | %d | %d / %d / %d | %d | %.0f | %s ± %s | %s ± %s |" % (s, st0[0].split()[-1] if st0 else "?", min(d["depth"] for d in ch),
                   len(glob.glob(f"{H}/runs/{s}/s*/est.out")), len(rows), sum(d["out"] == "D" for d in ch), sum(d["out"] == "C" for d in ch), sum(d["out"] == "B" for d in ch),
                   max(r[3] for r in rows), st.mean(d["fcpu"] + d["ccpu"] for d in ch), f"{st.mean(r[1] for r in rows):,.0f}", f"{se(1):,.0f}", f"{st.mean(r[2] for r in rows):,.0f}", f"{se(2):,.0f}"))
    return "\n".join(out)
t = rd("README.tmpl.md")
rep = {"ANSWER": rd("answer.md"), "NOTDONE": rd("notdone.md"),
       "D1LONG": nodetable(sorted(glob.glob(H + "/runs/L1_*/dnn.txt"))), "D1SHORT": nodetable(sorted(glob.glob(H + "/runs/s1_*/dnn.txt"))),
       "SPINE": nodetable(sorted(glob.glob(H + "/runs/L1_0/dnn.txt"))[:1]).split("\n")[0] + "\n" + "\n".join(
           x for x in nodetable([H + "/runs/L1_0/dnn.txt"] + sorted(glob.glob(H + "/runs/spine_*/dnn.txt"))).split("\n")[1:] if "| yes |" in x or x.startswith("|---")),
       "TWIN": nodetable(sorted(glob.glob(H + "/runs/twin32_*/dnn.txt"))), "DIVES": divetable(),
       "ESTIMATE": rd("estimate_text.md") + "\n\n```\n" + rd("logs/estimate.txt") + "\n```",
       "VALIDITY": "```\n" + rd("logs/validity.txt") + "\n```"}
for k, v in rep.items(): t = t.replace("@" + k + "@", v)
open(os.path.join(H, "README.md"), "w").write(t + "\n"); print("README.md", len(t.splitlines()), "lines")
