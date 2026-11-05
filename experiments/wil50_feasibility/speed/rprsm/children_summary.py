#!/usr/bin/env python3
"""Per-child table of run_children.sh: status, iterations, CPU and the
certified bound at iterations 50/100/200, base against the variant.
   usage: children_summary.py base name"""
import glob, os, re, sys
D = os.path.dirname(os.path.abspath(__file__))


def parse(path):
    """[(spec, result dict, trace)] for the parent (spec '-') and each child."""
    out, spec, tr = [], "-", {}
    for line in open(path):
        m = re.match(r"it\s+(\d+) lb (\S+)", line)
        if m: tr[int(m.group(1))] = float(m.group(2))
        if line.startswith("CHILD-SPEC"):
            spec = line.split()[1]
        if line.startswith("RESULT"):
            r = dict(kv.split("=", 1) for kv in line.split()[1:])
            out.append((spec, r, tr)); tr = {}
    return out


a, b = sys.argv[1], sys.argv[2]
tot = {a: [0, 0.0, 0, 0.0], b: [0, 0.0, 0, 0.0]}
print("%-6s %-6s  %-38s  %-38s" % ("parent", "child", a, b))
for f in sorted(glob.glob(os.path.join(D, "runs_child", a + "_p*.log"))):
    p = re.search(r"_p(\d+)\.log", f).group(1)
    g = f.replace("/%s_p" % a, "/%s_p" % b)
    if not os.path.exists(g): continue
    for (s1, r1, t1), (s2, r2, t2) in zip(parse(f), parse(g)):
        cells = []
        for name, r, t in ((a, r1, t1), (b, r2, t2)):
            st = "F" if r["status"] == "fathomed" else "-"
            cells.append("%s %4s it %6.1f s  lb100 %8.1f lb200 %8.1f" % (
                st, r["iters"], float(r["cpu"]), t.get(100, float("nan")), t.get(200, float("nan"))))
            if s1 != "-":
                k = 0 if st == "F" else 2
                tot[name][k] += int(r["iters"]); tot[name][k + 1] += float(r["cpu"])
        print("%-6s %-6s  %s  %s" % (p, s1, cells[0], cells[1]))
for name in (a, b):
    t = tot[name]
    print("%s: children fathomed %d it %.0f s; missed %d it %.0f s" % (name, t[0], t[1], t[2], t[3]))
