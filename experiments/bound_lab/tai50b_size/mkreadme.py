"""README.md from README.tmpl.md: the @NAME@ marks are filled from logs/ (run ./mkreport.sh first) and from the
hand-written answer.md and notdone.md.   python3 mkreadme.py"""
import os, glob, sys, statistics as st, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N, dives as DV, estimate as E
H = os.path.dirname(os.path.abspath(__file__)); rd = lambda f: open(os.path.join(H, f)).read().rstrip("\n")
def nodetable(files, bydepth=False):
    out = ["| node | on path | first solve: iterations | short of the cutoff | phase iterations | short after the phase | phase ended | outcome | CPU-s |", "|---|---|---|---|---|---|---|---|---|"]
    for d in sorted(N.load(files), key=lambda d: d["depth"]) if bydepth else N.load(files):
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
def stops():
    out = ["| set | seed | planned dives | finished | ended |", "|---|---|---|---|---|"]
    for d in sorted(glob.glob(H + "/runs/L[AB]*/s*")):
        if "LA0" in d: continue
        plan = open(d + "/estimate.par").read().split()[1]; c = open(d + "/cmd.txt").read()
        fin = len(open(d + "/dives.txt").readlines()) if os.path.exists(d + "/dives.txt") else 0
        end = "its count" if "exit 0" in c else ("its wall limit" if "exit 124" in c else "stopped by hand")
        out.append("| `%s` | %s | %s | %d | %s |" % (os.path.basename(os.path.dirname(d)), os.path.basename(d)[1:], plan if plan != "10000" else "by the clock", fin, end))
    return "\n".join(out)
def nbranch():
    n = b = 0
    for s_ in glob.glob(H + "/runs/L[AB]*"):
        if "LA0" in s_: continue
        for f, c in DV.dives(s_):
            if c[-1]["fath"]: n += 1; b += any(x["out"] == "B" for x in c)
    return n, b
t = rd("README.tmpl.md")
rep = {"ANSWER": rd("answer.md"), "NOTDONE": rd("notdone.md"),
       "D1LONG": nodetable(sorted(glob.glob(H + "/runs/L1_*/dnn.txt"))), "D1SHORT": nodetable(sorted(glob.glob(H + "/runs/s1_*/dnn.txt"))),
       "SPINE": nodetable(sorted(glob.glob(H + "/runs/L1_0/dnn.txt"))[:1]).split("\n")[0] + "\n" + "\n".join(
           x for x in nodetable([H + "/runs/L1_0/dnn.txt"] + sorted(glob.glob(H + "/runs/spine_*/dnn.txt")), True).split("\n")[1:] if "| yes |" in x or x.startswith("|---")),
       "TWIN": nodetable(sorted(glob.glob(H + "/runs/twin32_*/dnn.txt"))), "DIVES": divetable(),
       "ESTIMATE": rd("estimate_text.md") + "\n\n```\n" + rd("logs/estimate.txt") + "\n```",
       "VALIDITY": "```\n" + rd("logs/validity.txt") + "\n```"}
wc = rd("logs/webcheck.txt").split("\n")
tab = [x for x in wc if x.startswith("|")]
hi = [x.split("|") for x in tab[2:] if x.split("|")[1].strip().split()[0] in ("2", "3", "5")]
rep["WEBCHECK"] = "\n".join("  " + x for x in tab); rep["NCHECK"] = wc[0].split()[0]
rep["ABOVE2"] = str(sum(int(x[2]) for x in hi)); rep["WEBCOST"] = "```\n" + rd("logs/webcost.txt") + "\n```"
rep["STOPS"] = stops(); nd, nb = nbranch(); rep["NDIVES"] = str(nd); rep["NB1"] = str(nb)
for _ in range(2):
    for k, v in rep.items(): t = t.replace("@" + k + "@", v)
open(os.path.join(H, "README.md"), "w").write(t + "\n"); print("README.md", len(t.splitlines()), "lines")
