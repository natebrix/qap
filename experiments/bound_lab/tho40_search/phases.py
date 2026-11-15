"""Every cut phase in the cut-arm runs: depth, first-solve shortfall, pooled or not, result, iterations.
  python3 phases.py [ARM ...] [-q] [--path]    (default: every arm but base; the path-check runs only with --path)
Prints the list sorted by shortfall and a table by shortfall bin."""
import os, sys, glob, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, H
arms = [a for a in sys.argv[1:] if not a.startswith("-")]
rows = []
for d in sorted(glob.glob(os.path.join(H, "runs", "*_*"))):
    name = os.path.basename(d); arm = name.split("_", 1)[1]
    if arm == "base" or (arms and arm not in arms): continue
    if name.startswith("path") and "--path" not in sys.argv: continue   # incumbent 240,518 on the path: cannot fathom
    _, dnn = load(d)
    for k, e in enumerate(dnn):
        if e["cut"]:
            rows.append((e["cutoff"] - e["lb"], e["depth"], e["kept"] > 0, e["cstatus"] == "fathomed", e["cit"], e["cutoff"] - e["cbest"], e["iters"], name, k == 0, e["ccpu"]))
rows.sort()
if "-q" not in sys.argv:
    for r in rows:
        print("short %7.1f d%d %s -> %-4s after %4d it (%.0f s), left %7.1f; first solve %3d it  %s%s" % (r[0], r[1], "pool" if r[2] else "    ", "FATH" if r[3] else "miss", r[4], r[9], max(r[5], 0), r[6], r[7], (" ROOT(warm)" if r[7].startswith("w") else " ROOT(cold)") if r[8] else ""))
bins = [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3500, 5000, 1e9]
print("\nshortfall bin: phases, fathomed, median it (fathomed), median it (missed), median gain of misses")
def med(v): v = sorted(v); return v[len(v) // 2] if v else float("nan")
for lo, hi in zip(bins, bins[1:]):
    b = [r for r in rows if lo <= r[0] < hi]
    if not b: continue
    f = [r for r in b if r[3]]; m = [r for r in b if not r[3]]
    print("%5.0f-%-6.0f %3d %3d  %5.0f %5.0f  %7.1f" % (lo, min(hi, 99999), len(b), len(f), med([r[4] for r in f]), med([r[4] for r in m]), med([r[0] - r[5] for r in m])))
