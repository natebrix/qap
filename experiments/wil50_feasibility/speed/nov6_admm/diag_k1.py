#!/usr/bin/env python3
"""A dnntest -v -k 1 log (a certificate at every iteration): for the parent and each
warm child, the certified value at the first iterations (the warm start's transient
against the parent's final value), the iteration at which the running best crosses
the cutoff, and how the crossing predictor of DNN_ADMM_CROSS would have done from
the 10-iteration checks alone (its predicted iteration against the true one)."""
import re, sys, math
solves, cur, cutoff = [], None, None
for line in open(sys.argv[1], errors="replace"):
    if line.startswith("node ") or line.startswith("CHILD-PASS"):
        cur = {"it": {}, "best": {}}
        m = re.search(r"cutoff=(\S+)", line)
        if m: cutoff = float(m.group(1))
    m = re.match(r"it\s+(\d+) lb (\S+) \(cert ([0-9.eE+-]+)", line)
    if m and cur is not None:
        cur["best"][int(m.group(1))] = float(m.group(2)); cur["it"][int(m.group(1))] = float(m.group(3))
    if line.startswith("RESULT") and cur is not None:
        cur["result"] = line.strip(); solves.append(cur); cur = None
parent = solves[0]
pfinal = max(parent["best"].values())
print("cutoff %.2f; parent final certificate %.2f after %d iterations" % (cutoff, pfinal, max(parent["it"])))
for k, s in enumerate(solves[1:], 1):
    c = s["it"]; b = s["best"]
    cross = min((i for i in sorted(b) if b[i] > cutoff), default=None)
    first = " ".join("%d:%.0f" % (i, c[i]) for i in sorted(c) if i <= 12)
    print("child %d: cert at it 1-12: %s" % (k, first))
    print("         relative to parent's final: it1 %+.0f it5 %+.0f it10 %+.0f it20 %+.0f; crosses at %s" % (
        c[1] - pfinal, c.get(5, float("nan")) - pfinal, c.get(10, float("nan")) - pfinal,
        c.get(20, float("nan")) - pfinal, cross))
    if cross is None: continue
    # the predictor from the 10-iteration checks: at check it (multiple of 10) with the
    # running best below the cutoff, g = gain/it over the last interval, decayed by the
    # ratio of the last two gains; predicted crossing it + ceil((cutoff - b)/g)
    checks = [i for i in sorted(b) if i % 10 == 0 and i < cross]
    line_ = []
    for j, i in enumerate(checks):
        if j < 1: continue
        g = (b[i] - b[checks[j - 1]]) / 10.0
        if g <= 0: continue
        if j >= 2:
            g0 = (b[checks[j - 1]] - b[checks[j - 2]]) / 10.0
            if g0 > 0: g *= min(1.0, max(0.3, g / g0))
        need = (cutoff - b[i]) / g
        if need < 10: line_.append("at %d: predicted %d (true %d)" % (i, i + math.ceil(need), cross))
    print("         predictor: %s; gains per it over the last 20 before the crossing: %s" % (
        "; ".join(line_) if line_ else "no prediction before the regular check at %d" % (10 * math.ceil(cross / 10.0)),
        " ".join("%.1f" % (b[i] - b[i - 1]) for i in range(max(2, cross - 19), cross + 1) if i in b and i - 1 in b)))
