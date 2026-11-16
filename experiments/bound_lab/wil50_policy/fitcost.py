"""The hand-over cost on wil50/sko49: what the children of a phase that missed cost, by the gap it left.
  python3 fitcost.py N UNITCUTOFF trace.txt [trace.txt ...]
Prints each missed phase (m, gap left, descendants, iterations) binned by gap, and a least-squares fit of
iterations / m = N0 + A g + B g^2, g in 0.1% of the cutoff."""
import sys, numpy as np
from tracelib import load, below
N = int(sys.argv[1]); U = float(sys.argv[2]) * 0.001
rows = []
for f in sys.argv[3:]:
    nodes = load(f); b = below(nodes)
    for e, (k, it) in zip(nodes, b):
        if e["cut"] and not e["fathomed"]:
            rows.append((N - e["depth"], e["left"], k, it, e["cit"], e["cstatus"], e["depth"], e["short"]))
rows.sort(key=lambda r: r[1])
g = np.array([r[1] / U for r in rows]); y = np.array([r[3] / r[0] for r in rows])
print("missed phases: %d" % len(rows))
edges = [0, 5, 10, 20, 30, 50, 75, 100, 150, 250, 1e9]
print("gap left (units) | phases | phase it mean | children it per child: mean, median | DNN nodes below per child")
for lo, hi in zip(edges, edges[1:]):
    s = [r for r in rows if lo <= r[1] < hi]
    if not s: continue
    pc = [r[3] / r[0] for r in s]
    print("%4.0f-%-4.0f | %3d | %5.0f | %7.1f %7.1f | %5.2f" % (lo, min(hi, 999), len(s), np.mean([r[4] for r in s]), np.mean(pc), np.median(pc), np.mean([r[2] / r[0] for r in s])))
X = np.vstack([np.ones_like(g), g, g * g]).T
for name, sel in [("all", g < 1e9), ("gap <= 150 units", g * U <= 150)]:
    c, *_ = np.linalg.lstsq(X[sel], y[sel], rcond=None)
    # weighted in log space is more honest for the scatter; report plain LS and the residual factor
    pred = X[sel] @ c
    print("fit (%s): N0 %.1f A %.2f B %.3f | rms log10 residual %.2f" % (name, c[0], c[1], c[2], np.sqrt(np.mean((np.log10(np.maximum(pred, 1)) - np.log10(np.maximum(y[sel], 1))) ** 2))))
