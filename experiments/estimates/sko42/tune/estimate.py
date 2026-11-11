#!/usr/bin/env python3
"""sko42's tree estimate from stage B's sample and complete subtrees ("anchors"),
in DNN nodes and in CPU:

  estimate.py [--table sampling/stageB_summary.md] [--top 451] [--top-cpu 40762]
              [--front 10726] ANCHOR_DIR=SHORTFALL ...

Each ANCHOR_DIR is a subtree.sh output (its summary.txt: "N DNN nodes ... DNN CPU X s")
whose node's give-up shortfall in stage B is SHORTFALL. Unlike wil50's refit.py,
subtrees are indexed by the give-up shortfall (sko42's converge step moved the
farthest four by 3-5 units) and their CPU is measured, not nodes times a fixed price:
a node that gets a cut phase costs 100-500 s, one that does not 5-45 s.

For each sampled frontier node, s its give-up shortfall:
  s <= 0 (fathomed in the first solve)  1 node, its own first-solve CPU (the table's)
  s > 0                                  the anchors' size at s, under two models:
    interp   log-linear between the anchors by s, flat below the smallest, the fit's
             slope above the largest
    fit      log(y) = a + b s by least squares over the anchors, times exp(sigma^2/2)
             (the mean when sizes scatter log-normally about the fit)
Totals: the --top nodes and --top-cpu CPU above the frontier (stage A's) plus
--front frontier nodes times the sample mean; a 95% bootstrap over the sample."""
import math, os, random, re, sys

H = os.path.dirname(os.path.abspath(__file__))


def sample(table):
    rows = []
    for ln in open(table):
        if not ln.startswith("| ") or ln.startswith("| node") or ln.startswith("|---"):
            continue
        c = [x.strip() for x in ln.strip().strip("|").split("|")]
        if len(c) < 15:
            continue
        try:
            short, cpu = float(c[4]), float(c[7])
        except ValueError:
            continue
        rows.append((-1.0 if c[5] == "fathomed" else short, cpu))
    return rows


def anchor(d):
    t = open(os.path.join(d, "summary.txt")).readline()
    n = float(re.search(r"([\d,]+) DNN nodes", t).group(1).replace(",", ""))
    cpu = float(re.search(r"DNN CPU ([\d,.]+) s", t).group(1).replace(",", ""))
    return n, cpu


def fit(pts):
    xs = [x for x, _ in pts]; ys = [math.log(y) for _, y in pts]; k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    s2 = sum((y - a - b * x) ** 2 for x, y in zip(xs, ys)) / max(1, k - 2)
    return a, b, s2


def interp(pts, b):
    p = sorted(pts)
    def f(s):
        if s <= p[0][0]:
            return p[0][1]
        for (s0, y0), (s1, y1) in zip(p, p[1:]):
            if s <= s1:
                return math.exp(math.log(y0) + (math.log(y1) - math.log(y0)) * (s - s0) / (s1 - s0))
        return p[-1][1] * math.exp(b * (s - p[-1][0]))
    return f


def main():
    args = sys.argv[1:]
    table, top, top_cpu, front = os.path.join(H, "sampling", "stageB_summary.md"), 451, 40762.0, 10726
    while args and args[0].startswith("--"):
        o = args.pop(0)
        if o == "--table": table = args.pop(0)
        elif o == "--top": top = int(args.pop(0))
        elif o == "--top-cpu": top_cpu = float(args.pop(0))
        elif o == "--front": front = int(args.pop(0))
        else: raise SystemExit(f"unknown option {o}")
    anchors = []
    for a in args:
        d, _, s = a.rpartition("=")
        n, cpu = anchor(d)
        anchors.append((float(s), n, cpu, os.path.basename(d.rstrip("/"))))
    if len(anchors) < 2:
        raise SystemExit("two anchors or more")
    S = sample(table)
    miss = [s for s, _ in S if s > 0]
    print(f"stage B: {len(S)} sampled, {len(S) - len(miss)} fathomed in the first solve, misses {min(miss):.1f}-{max(miss):.1f} short;"
          f" frontier {front:,}, {top:,} nodes and {top_cpu:,.0f} CPU s above it")
    print("anchors (give-up short, DNN nodes, DNN CPU s): " + "; ".join(f"{l} {s:.1f} {n:,.0f} {c:,.0f}" for s, n, c, l in sorted(anchors)))
    rng = random.Random(1)
    print("| quantity | model | fit: per 100 short, residual factor | mean per frontier node | total (95%) |")
    print("|---|---|---|---|---|")
    for qi, q in ((1, "DNN nodes"), (2, "CPU s")):
        pts = [(a[0], a[qi]) for a in anchors]
        a, b, s2 = fit(pts)
        models = [("interp", interp(pts, b)), ("fit", lambda s: math.exp(a + b * s + s2 / 2))]
        for name, f in models:
            z = [(1.0 if qi == 1 else c) if s <= 0 else f(s) for s, c in S]
            m = sum(z) / len(z)
            boots = sorted(sum(z[rng.randrange(len(z))] for _ in z) / len(z) for _ in range(4000))
            base = top if qi == 1 else top_cpu
            tot = base + front * m
            lo, hi = base + front * boots[100], base + front * boots[3899]
            tag = f"x{math.exp(100 * b):.1f}, x{math.exp(math.sqrt(s2)):.1f}" if name == "fit" else ""
            unit = "" if qi == 1 else f" = {tot / 3600:,.0f} core-hours"
            print(f"| {q} | {name} | {tag} | {m:,.1f} | {tot:,.0f} ({lo:,.0f}-{hi:,.0f}){unit} |")


if __name__ == "__main__":
    main()
