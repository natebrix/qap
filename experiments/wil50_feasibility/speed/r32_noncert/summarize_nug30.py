#!/usr/bin/env python3
"""Summarize nug30/ of run_nug30.sh (base_*, nc_*) and the second pass in
reverse order (base2_*, nc2_*): outcome, iterations, ms per iteration, the
certified bound at the last check, the largest certificate margin, form R."""
import glob, os, re, statistics
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nug30")


def parse(f):
    r = {"margins": []}
    for line in open(f):
        m = re.match(r"it\s+(\d+) lb \S+ \(cert \S+, margin (\S+?)[,)]", line)
        if m: r["margins"].append(float(m.group(2)))
        m = re.match(r"form R \(ms per iteration\): F (\S+)\s+syrk (\S+)\s+symmetrize (\S+)", line)
        if m: r["formR"] = sum(float(x) for x in m.groups())
        if line.startswith("RESULT"):
            for kv in line.split()[1:]:
                k, v = kv.split("="); r[k] = v
    return r if "status" in r else None


tags = sorted({re.match(r"base_(.+)\.log", os.path.basename(f)).group(1)
               for f in glob.glob(os.path.join(D, "base_*.log"))})
tot = {c: [0.0, 0] for c in ("base", "nc", "base2", "nc2")}
print("%-12s %-22s %-22s %8s %8s %9s %9s" % ("node", "base: status it ms/it", "nc: status it ms/it",
                                           "lb diff", "pass2", "formR b", "formR nc"))
for t in tags:
    r = {c: parse(os.path.join(D, "%s_%s.log" % (c, t))) for c in tot}
    for c in tot:
        if r[c]: tot[c][0] += float(r[c]["cpu"]); tot[c][1] += int(r[c]["iters"])
    def s(x): return "%s %4s %5.1f" % (x["status"][:4], x["iters"], float(x["ms_per_iter"]))
    p2 = (float(r["nc2"]["ms_per_iter"]) / float(r["base2"]["ms_per_iter"])) if r["nc2"] and r["base2"] else float("nan")
    print("%-12s %-22s %-22s %+8.4f %8.4f %9.2f %9.2f" % (
        re.sub(r"(\d+)_(\d+)_?", r"\1:\2,", t).rstrip(","), s(r["base"]), s(r["nc"]),
        float(r["nc"]["lb"]) - float(r["base"]["lb"]), p2, r["base"].get("formR", 0), r["nc"].get("formR", 0)))
for c in tot:
    if tot[c][1]: print("%-6s %6d iterations, %7.1f s, %.2f ms per iteration" % (c, tot[c][1], tot[c][0], tot[c][0] / tot[c][1] * 1e3))
mb = [x for f in glob.glob(os.path.join(D, "base_*.log")) for x in parse(f)["margins"]]
mn = [x for f in glob.glob(os.path.join(D, "nc_*.log")) for x in parse(f)["margins"]]
print("certificate margins: base max %.2e median %.2e; nc max %.2e median %.2e" % (
    max(mb), statistics.median(mb), max(mn), statistics.median(mn)))
