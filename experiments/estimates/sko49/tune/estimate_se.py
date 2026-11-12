#!/usr/bin/env python3
"""sko49's tree estimate with the anchors' own error in the interval:

  estimate_se.py [estimate.py's options] ANCHOR_DIR=SHORTFALL ...

`../../sko42/tune/estimate.py` reads the anchors as exact and bootstraps only
stage B's sample. On sko42 the anchors were complete subtrees; here most are
probes (`../../probe.py`), each with a standard error of up to 90% of its value,
and the tail they cover carries most of the total. This takes the same inputs
and prints, for CPU s and DNN nodes:
  interp  estimate.py's interp model. The 95% interval draws each probe's value
          from a log-normal with the probe's mean and SE (a complete subtree is
          exact), and stage B's sample with replacement, jointly.
  iso     interp through the anchors' isotonic fit: anchors sorted by shortfall,
          and any anchor dearer than the next pooled with it (pool adjacent
          violators; a block is its anchors' mean shortfall and mean cost).
          Replicates at nearly the same shortfall (Cp326 and Cr338: 349.1 and
          350.2 short, 170,185 and 27,433 CPU s) make interp swing; this pools
          them on the prior that the mean cost does not fall as the shortfall
          grows. Beyond the last block it extrapolates with the fit's slope, as
          interp does. The interval resamples the anchors as fit does.
  fit     estimate.py's fit model (log y = a + b s, times exp(sigma^2/2)). The
          interval resamples the anchors with replacement (the fit's scatter
          already holds the probes' error) and stage B's sample, jointly.
The interp and fit point estimates equal estimate.py's."""
import math, os, random, re, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "sko42", "tune"))
import estimate as E  # noqa: E402


def anchor_se(d):
    """(nodes, nodes SE, CPU, CPU SE) from a summary.txt; SE 0 when the line has none."""
    t = open(os.path.join(d, "summary.txt")).readline()
    n, cpu = E.anchor(d)
    m = re.search(r"DNN CPU [\d,.]+ s \(SE ([\d,.]+) s; nodes SE ([\d,.]+)\)", t)
    cse, nse = (float(m.group(1).replace(",", "")), float(m.group(2).replace(",", ""))) if m else (0.0, 0.0)
    return n, nse, cpu, cse


def lognormal(rng, mean, se):
    if se <= 0 or mean <= 0:
        return mean
    s2 = math.log(1 + (se / mean) ** 2)
    return math.exp(rng.gauss(math.log(mean) - s2 / 2, math.sqrt(s2)))


def total(S, f, qi, base, front):
    z = [(1.0 if qi == 1 else c) if s <= 0 else f(s) for s, c in S]
    return base + front * sum(z) / len(z)


def pava(pts):
    """the isotonic (non-decreasing) fit of (s, y) as (mean s, mean y) blocks"""
    blocks = []                                     # [sum s, sum y, count]
    for s, y in sorted(pts):
        blocks.append([s, y, 1])
        while len(blocks) > 1 and (blocks[-2][1] / blocks[-2][2] >= blocks[-1][1] / blocks[-1][2]
                                   or blocks[-2][0] / blocks[-2][2] >= blocks[-1][0] / blocks[-1][2]):
            s1, y1, n1 = blocks.pop()
            blocks[-1][0] += s1; blocks[-1][1] += y1; blocks[-1][2] += n1
    return [(bs / n, by / n) for bs, by, n in blocks]


def models(pts):
    """(interp, iso, fit) as functions of the shortfall"""
    a, b, s2 = E.fit(pts)
    return E.interp(pts, b), E.interp(pava(pts), b), (lambda s: math.exp(a + b * s + s2 / 2))


def main():
    args = sys.argv[1:]
    table, top, top_cpu, front = None, None, None, None
    while args and args[0].startswith("--"):
        o = args.pop(0)
        if o == "--table": table = args.pop(0)
        elif o == "--top": top = int(args.pop(0))
        elif o == "--top-cpu": top_cpu = float(args.pop(0))
        elif o == "--front": front = int(args.pop(0))
        else: raise SystemExit(f"unknown option {o}")
    if None in (table, top, top_cpu, front):
        raise SystemExit("--table, --top, --top-cpu and --front are all needed")
    A = []
    for a in args:
        d, _, s = a.rpartition("=")
        A.append((float(s), *anchor_se(d), os.path.basename(d.rstrip("/"))))
    if len(A) < 3:
        raise SystemExit("three anchors or more")
    S = E.sample(table)
    rng = random.Random(1)
    B = 2000
    print(f"{len(A)} anchors, {len(S)} sampled frontier nodes, {B} draws a model")
    print("| quantity | model | total | 95% (sample only) | 95% (sample and anchors) | cloud core-hours (95%) |")
    print("|---|---|---|---|---|---|")
    for qi, q, base in ((1, "DNN nodes", top), (3, "CPU s", top_cpu)):
        pts = [(x[0], x[qi]) for x in A]
        if qi == 3:
            print("iso blocks (CPU s): " + ", ".join(f"{bs:.1f}: {by:,.0f}" for bs, by in pava(pts)))
        for mi, (name, f) in enumerate(zip(("interp", "iso", "fit"), models(pts))):
            point = total(S, f, qi, base, front)
            only, both = [], []
            for _ in range(B):
                Sb = [S[rng.randrange(len(S))] for _ in S]
                only.append(total(Sb, f, qi, base, front))
                if name == "interp":
                    pb = [(x[0], lognormal(rng, x[qi], x[qi + 1])) for x in A]
                    g = models(pb)[0]
                else:
                    while True:
                        ab = [A[rng.randrange(len(A))] for _ in A]
                        if len({x[0] for x in ab}) >= 3:
                            break
                    g = models([(x[0], x[qi]) for x in ab])[mi]
                both.append(total(Sb, g, qi, base, front))
            only.sort(); both.sort()
            lo, hi = int(0.025 * B), int(0.975 * B) - 1
            ch = f"{point / 3600:,.0f} ({both[lo] / 3600:,.0f}-{both[hi] / 3600:,.0f})" if qi == 3 else ""
            print(f"| {q} | {name} | {point:.2e} | {only[lo]:.2e}-{only[hi]:.2e} | {both[lo]:.2e}-{both[hi]:.2e} | {ch} |")


if __name__ == "__main__":
    main()
