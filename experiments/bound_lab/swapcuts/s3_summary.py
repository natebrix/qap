#!/usr/bin/env python3
"""Stage 3 tables from the joint-rounds logs: logs/s3_<node>_jT*.log (triangles alone) against
logs/s3_<node>_jTX*.log (triangles + exchange cuts), certified bound after each round.
  python3 s3_summary.py > logs/s3_summary.txt
"""
import glob, os, re
here = os.path.dirname(os.path.abspath(__file__))


def rounds(path):
    out = {}; fin = None; cmd = ""
    for line in open(path):
        if line.startswith("# "):
            cmd = line[2:].strip()
        m = re.match(r"round +(\d+) +lb (-?[\d.]+) .*status (\w+) +iters (\d+)", line)
        if m:
            out[int(m.group(1))] = (float(m.group(2)), m.group(3), int(m.group(4)))
        m = re.match(r"quiet round +lb (-?[\d.]+) .*status (\w+) +iters (\d+)", line)
        if m:
            out["quiet"] = (float(m.group(1)), m.group(2), int(m.group(3)))
        if line.startswith("FINAL"):
            fin = line.strip()
    return out, fin, cmd


def viol(path):
    v = []
    for line in open(path):
        m = re.search(r"exchange separation (\d+): (\d+) of (\d+) cuts violated .*sum of violations ([\d.]+)", line)
        if m:
            v.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4))))
    return v


for t in sorted(glob.glob(os.path.join(here, "logs", "s3_*_jT*.log"))):
    name = os.path.basename(t)
    if "_jTX" in name:
        continue
    x = t.replace("_jT", "_jTX")
    if not os.path.exists(x):
        continue
    rt, ft, ct = rounds(t)
    rx, fx, cx = rounds(x)
    print("== %s  vs  %s" % (name, os.path.basename(x)))
    print("   T : " + ct)
    print("   TX: " + cx)
    cutoff = None
    m = re.search(r"-u (\d+)", ct)
    if m:
        cutoff = float(m.group(1))
    its = 0
    print("   %-6s %10s %12s %12s %8s   %s" % ("round", "iters", "T", "TX", "TX - T", "(shortfall of T / TX to the best known value)" if cutoff else ""))
    keys = [k for k in rt if k != "quiet"]
    for k in sorted(set(keys) | set(k for k in rx if k != "quiet")):
        a = rt.get(k); b = rx.get(k)
        its += (a or b)[2]
        print("   %-6s %10d %12s %12s %8s   %s" % (k, its, "%.1f" % a[0] if a else "-", ("%.1f" % b[0] if b else "-") + ("*" if b and b[1] == "fathomed" else ""),
              "%+.1f" % (b[0] - a[0]) if a and b else "", "%.0f / %.0f" % (cutoff - a[0], cutoff - b[0]) if cutoff and a and b else ""))
    if "quiet" in rt or "quiet" in rx:
        a = rt.get("quiet"); b = rx.get("quiet")
        print("   %-6s %10s %12s %12s %8s" % ("quiet", "+%d" % (a or b)[2], "%.1f" % a[0] if a else "-", "%.1f" % b[0] if b else "-", "%+.1f" % (b[0] - a[0]) if a and b else ""))
    best = lambda r: max(v[0] for v in r.values())
    if rt and rx:
        print("   best certified: T %.1f, TX %.1f, TX - T %+.1f%s" % (best(rt), best(rx), best(rx) - best(rt),
              "   (* = fathomed, the run stops there)" if any(v[1] == "fathomed" for v in rx.values()) else ""))
    v = viol(x)
    if v:
        print("   exchange cuts violated at TX's iterates (round: count, sum of violations): " +
              "  ".join("%d: %d, %.0f" % (r, c, s) for r, c, n, s in v))
    print()
