#!/usr/bin/env python3
"""Part 2 table from logs/p2pl_*.log (plain DNN, dnntest with tho40.par), logs/p2cp_*.log (CP phase
with the cutoff at 241,500: the fathom iteration is the first check whose certified value passes
240,514.05) and logs/p2u_*.log (CP phase with -u 240516, where run)."""
import os, re, sys, math
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CUT = 240514.05
def rd(tag):
    p = os.path.join(D, "logs", tag + ".log")
    return open(p).read() if os.path.exists(p) else None
nodes = [l.rstrip("\n").split("\t") for l in open(os.path.join(D, "nodes.tsv"))][1:]
out = []
for tag, d, m, spec, pdnn, pshort, kids, sub, subcpu in nodes:
    r = dict(tag=tag, d=int(d), m=int(m), spec=spec, pshort=pshort, sub=sub)
    t = rd("p2pl_" + tag)
    if t:
        its = [(int(i), float(v)) for i, v in re.findall(r"^it +(\d+) lb (\S+)", t, re.M)]
        best = -1e30; tr = []
        for i, v in its:
            best = max(best, v); tr.append((i, best))
        at = lambda n: max([b for i, b in tr if i <= n] or [float("nan")])
        r["pl300"] = at(300); r["pl600"] = at(600)
        res = re.search(r"RESULT .* lb=(\S+) .*status=(\S+) iters=(\d+)", t)
        if res:
            r["plfin"] = float(res.group(1)); r["plst"] = res.group(2); r["plit"] = int(res.group(3))
        elif tr:
            r["plfin"] = tr[-1][1]; r["plst"] = "running"; r["plit"] = tr[-1][0]
        # rate over the last 300 iterations
        if tr:
            last = tr[-1][0]; r["pllast300"] = tr[-1][1] - at(last - 300)
    t = rd("p2cp_" + tag)
    if t:
        r0 = re.search(r"^round  0  lb (\S+)", t, re.M)
        if r0: r["c300"] = float(r0.group(1))
        certs = [(int(i), float(c)) for i, c in re.findall(r"^cp it +(\d+) cert (\S+) ", t, re.M)]
        if certs:
            best = -1e30; tr = []
            for i, v in certs:
                best = max(best, v); tr.append((i, best))
            f = [i for i, b in tr if b > CUT]
            r["cpfath"] = f[0] if f else None
            r["cpfin"] = tr[-1][1]; r["cpit"] = tr[-1][0]
            r["cplast200"] = tr[-1][1] - max([b for i, b in tr if i <= tr[-1][0] - 200] or [float("nan")])
            for n in (100, 200, 400, 800):
                v = [b for i, b in tr if i <= n]
                r["cp%d" % n] = v[-1] if v and tr[-1][0] >= n else float("nan")
            pool = re.findall(r"\| pool (\d+) sep", t)
            r["pool"] = int(pool[-1]) if pool else 0
    t = rd("p2u_" + tag)
    if t:
        fin = re.search(r"^cp solve  lb (\S+) .*status (\S+)  iters (\d+)", t, re.M)
        if fin: r["u"] = (float(fin.group(1)), fin.group(2), int(fin.group(3)))
    out.append(r)
g = lambda r, k, f="%.1f": (f % r[k]) if k in r and r[k] is not None and not (isinstance(r[k], float) and math.isnan(r[k])) else "-"
print("| node | depth | m | production shortfall (subtree nodes) | plain@300 short | plain end short (its, last-300 gain) | CP: fathom at | CP value at 100/200/400/800 | CP end short (gain/200 its) | CP gain vs plain end | share of plain end shortfall | -u 240516 run |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in out:
    pl300 = CUT - r["pl300"] if "pl300" in r else None
    plfin = CUT - r["plfin"] if "plfin" in r else None
    cpfin = CUT - r["cpfin"] if "cpfin" in r else None
    gain = (r["cpfin"] - r["plfin"]) if ("cpfin" in r and "plfin" in r) else None
    share = gain / plfin if gain is not None and plfin and plfin > 0 else None
    fath = ("it %d (300+%d)" % (r["cpfath"], r["cpfath"])) if r.get("cpfath") else ("no" if "cpfin" in r else "-")
    u = "%s %.2f @%d" % (r["u"][1], r["u"][0], r["u"][2]) if "u" in r else "-"
    print("| %s `%s` | %d | %d | %s (%s) | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        r["tag"], r["spec"], r["d"], r["m"], r["pshort"], r["sub"],
        "%.0f" % pl300 if pl300 is not None else "-",
        ("%.0f (%s %d, %+.0f)" % (plfin, r["plst"], r["plit"], r.get("pllast300", float("nan")))) if plfin is not None else "-",
        fath,
        "/".join("%.0f" % (r["cp%d" % n] - CUT) if ("cp%d" % n) in r and not math.isnan(r["cp%d" % n]) else "-" for n in (100, 200, 400, 800)) if "cpfin" in r else "-",
        ("%.0f (@%d, %+.0f)" % (cpfin, r["cpit"], r["cplast200"])) if cpfin is not None else "-",
        "%+.0f" % gain if gain is not None else "-",
        ("%.2f" % share) if share is not None else ("plain fathoms" if plfin is not None and plfin <= 0 else "-"), u))
import json
json.dump(out, open(os.path.join(D, "part2.json"), "w"), indent=1)
