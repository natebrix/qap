"""ms per DNN iteration by arm and depth: unpooled first solves, pooled first solves, phases; and the
pooled children of unfathomed phases against the same nodes in base.   python3 costs.py ARM..."""
import sys, collections
from summarize import load
key = lambda e: frozenset(e["fixes"].split())
_, b = load("base"); bk = {key(e): e for e in b}
for arm in sys.argv[1:]:
    _, c = load(arm)
    t = collections.defaultdict(lambda: [0.0, 0] * 3)
    for e in c:
        x = t[e["depth"]]; j = 2 if e["kept"] else 0
        x[j] += e["cpu"]; x[j + 1] += e["iters"]; x[4] += e["ccpu"]; x[5] += e["cit"]
    print(f"\n{arm}: ms per iteration | depth | unpooled first solve | pooled first solve | phase |")
    for d in sorted(t):
        x = t[d]; f = lambda a, n: "%.1f (%s it)" % (1000 * a / n, f"{n:,}") if n else "-"
        print(f"| {d} | {f(x[0], x[1])} | {f(x[2], x[3])} | {f(x[4], x[5])} |")
    # pooled nodes that base also has: outcome and iterations
    for d in sorted(t):
        m = [(e, bk[key(e)]) for e in c if e["depth"] == d and e["kept"] and key(e) in bk]
        if not m: continue
        print("  depth %d pooled nodes also in base: %d; fathomed in the first solve %d (base %d); first-solve iterations %s (base %s); CPU s %.0f (base %.0f); median shortfall of the unfathomed %.1f (base %.1f)" % (
            d, len(m), sum(e["status"] == "fathomed" for e, _ in m), sum(g["status"] == "fathomed" for _, g in m),
            f'{sum(e["iters"] for e, _ in m):,}', f'{sum(g["iters"] for _, g in m):,}', sum(e["cpu"] for e, _ in m), sum(g["cpu"] for _, g in m),
            (lambda v: sorted(v)[len(v) // 2] if v else float("nan"))([e["cutoff"] - e["lb"] for e, _ in m if e["status"] != "fathomed"]),
            (lambda v: sorted(v)[len(v) // 2] if v else float("nan"))([g["cutoff"] - g["lb"] for _, g in m if g["status"] != "fathomed"])))
