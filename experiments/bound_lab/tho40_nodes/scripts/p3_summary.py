#!/usr/bin/env python3
"""Part 3 tables: the rounds form (logs/p3r_ARM_NODE.log) per round, and the CP form with
round 0's boxes in the pool (logs/p3cpb_NODE.log) against the triangles-only CP run (p2cp)."""
import glob, os, re
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")
CUT = 240514.05
nodes = sorted({os.path.basename(f)[:-4].split("_", 2)[2] for f in glob.glob(os.path.join(D, "p3r_*.log"))})
for nd in nodes:
    arms = {}
    for arm in ("tri", "five", "seven", "box"):
        p = os.path.join(D, "p3r_%s_%s.log" % (arm, nd))
        if not os.path.exists(p): continue
        t = open(p).read()
        rs = re.findall(r"^round +(\d+)  lb (\S+) .*?iters (\d+)  cpu (\S+)(.*)$", t, re.M)
        arms[arm] = rs
    if not arms: continue
    print("\nRounds form, node %s (shortfall to 240,514.05 after each round; iterations cumulative):\n" % nd)
    names = [a for a in ("tri", "five", "seven", "box") if a in arms]
    print("| round | iterations | " + " | ".join({"tri": "triangles", "five": "+ five 500", "seven": "+ seven 200", "box": "+ box 200"}[a] for a in names) + " |")
    print("|---|---|" + "---|" * len(names))
    nr = max(len(arms[a]) for a in names)
    for k in range(nr):
        its = None; cells = []
        for a in names:
            if k < len(arms[a]):
                cum = sum(int(x[2]) for x in arms[a][:k + 1]); its = its or cum
                best = max(float(x[1]) for x in arms[a][:k + 1])
                cells.append("%.0f" % (CUT - best))
            else: cells.append("-")
        print("| %d | %s | %s |" % (k, its, " | ".join(cells)))
    for a in names:
        last = arms[a][-1][4]
        ex = re.findall(r"(five|seven|box) (\d+) \(max ([0-9.]+)\) added\w+ (\d+)", last)
        tri = re.search(r"active (\d+)  violated (\d+) \(max ([0-9.]+)\)", last)
        cpu = sum(float(x[3]) for x in arms[a])
        print("- %s: last round active %s, triangles violated %s (max %s)%s; cpu %.0f s" % (a, tri.group(1) if tri else "?", tri.group(2) if tri else "?", tri.group(3) if tri else "?",
              "".join("; %s violated %s (max %s)" % (e[0], e[1], e[2]) for e in ex), cpu))
print("\nCP form: triangles (p2cp) against triangles + round 0's 200 boxes (p3cpb), certified shortfall to 240,514.05:\n")
print("| node | CP it 100 | 200 | 400 | 600 | 800 | pool at the end | at the give-up point: violated (max) | at the CP run's last iterate: violated (max) |")
print("|---|---|---|---|---|---|---|---|---|")
for f in sorted(glob.glob(os.path.join(D, "p3cpb_*.log"))):
    nd = os.path.basename(f)[6:-4]
    for tag, p in (("triangles", os.path.join(D, "p2cp_%s.log" % nd)), ("+ boxes", f)):
        if not os.path.exists(p): continue
        t = open(p).read()
        certs = {int(i): float(c) for i, c in re.findall(r"^cp it +(\d+) cert (\S+) ", t, re.M)}
        best = {}; b = -1e30
        for i in sorted(certs): b = max(b, certs[i]); best[i] = b
        pool = re.findall(r"\| pool (\d+) sep", t)
        r0 = re.findall(r"^round  0  lb .*$", t, re.M)
        def viol(line):
            tri = re.search(r"violated (\d+) \(max ([0-9.]+)\)", line)
            ex = re.findall(r"(five|seven|box) (\d+) \(max ([0-9.]+)\)", line)
            return "tri %s (%s)" % (tri.group(1), tri.group(2)) + "".join(", %s %s (%s)" % e for e in ex) if tri else "-"
        print("| %s %s | %s | %s | %s |" % (nd, tag, " | ".join("%.0f" % (CUT - best[n]) if n in best else "-" for n in (100, 200, 400, 600, 800)),
              pool[-1] if pool else "-", viol(r0[0]) if r0 else "-") + (" %s |" % (viol(r0[1]) if len(r0) > 1 else "-")))
