"""Replay stop rules offline on recorded phase histories (QAP_DNN_TRACE_HIST: running best every 10 it),
after tho40_search/stop/replay_rules.py.
  python3 replay.py INST RUN [RUN ...]     INST wil50 | sko49; RUN = runs/ directory names
Cost of a phase under a rule = iterations to its stop + (if it did not fathom) the children's cost at the gap
left, m x c(gap), c the per-child cost measured on missed phases whose children ran complete (fitcost.py;
the table CHILD below, interpolated in the log). A rule more patient than the recorded stop is censored: it
keeps the recorded end. m = N - depth."""
import sys, os, math
from tracelib import load
H = os.path.dirname(os.path.abspath(__file__))
INST = sys.argv[1]
N, CUT = (50, 48816.0) if INST == "wil50" else (49, 23386.0)
UNIT = 0.001 * CUT
# per-child iterations after a miss, by gap left in 0.1% of the cutoff: wil50 a251+a325+a220 laptop traces
# (364 missed phases, fitcost.py), bins' means at the bins' centres
CHILD = [(0.05, 32.0), (0.15, 34.5), (0.31, 40.3), (0.51, 60.4), (0.82, 131), (1.28, 301), (1.79, 697), (2.56, 1685), (4.1, 5607), (5.5, 22670)]
def child(g):
    if g <= CHILD[0][0]: return CHILD[0][1]
    for (g0, c0), (g1, c1) in zip(CHILD, CHILD[1:]):
        if g <= g1: return math.exp(math.log(c0) + (math.log(c1) - math.log(c0)) * (g - g0) / (g1 - g0))
    return CHILD[-1][1] * (g / CHILD[-1][0]) ** 4
def prev(h, it, W, g0):
    for jt, gj in h:
        if jt > it - W: break
        g0 = gj
    return g0
def stall(W, f=3.0):
    def r(p, h, k):
        it, g = h[k]
        return it >= W and g > f * max(0.0, prev(h, it, W, p["short"]) - g)
    return r
def giveup(pairs, base):
    return lambda p, h, k: any(h[k][0] >= T and h[k][1] > G for T, G in pairs) or base(p, h, k)
def cap(T, base):
    return lambda p, h, k: h[k][0] >= T or base(p, h, k)
def handover(n0, a, b, phi=1.5, W=200):
    def r(p, h, k):
        it, g = h[k]
        if it < W: return False
        rate = max(0.0, prev(h, it, W, p["short"]) - g) / W; u = g / UNIT
        return phi * g > rate * p["m"] * (n0 + a * u + b * u * u)
    return r
base = stall(600)
RULES = [("recorded: cut.stall 3 x cut.iters 200", base), ("stall W 400", stall(400)), ("stall W 300", stall(300)), ("stall W 200", stall(200)),
         ("cut.cp_iters 1500", cap(1500, base)), ("cut.cp_iters 1000", cap(1000, base)),
         ("cut.giveup 400:40,800:15", giveup([(400, 40 * UNIT / 48.8), (800, 15 * UNIT / 48.8)], base)),
         ("handover 32:0:250", handover(32, 0, 250)), ("handover 32:0:250 phi 1", handover(32, 0, 250, 1.0)), ("handover 32:0:250 phi 2.5", handover(32, 0, 250, 2.5)),
         ("handover 32:0:250 W 100", handover(32, 0, 250, 1.5, 100)), ("handover 32:0:250 W 300", handover(32, 0, 250, 1.5, 300)),
         ("handover 32:20:100", handover(32, 20, 100)), ("handover 32:0:60", handover(32, 0, 60)), ("handover 50:0:250", handover(50, 0, 250)),
         ("handover 20:2:3 (tho40's)", handover(20, 2, 3))]
ph = []
for run in sys.argv[2:]:
    for e in load(os.path.join(H, "runs", run, "dnn.txt")):
        if e["cut"] and e["hist"]:
            ph.append(dict(run=run, m=N - e["depth"], short=e["short"], hist=[(t, e["cutoff"] - b) for t, b in e["hist"]], its=e["cit"], fath=e["cstatus"] == "fathomed", left=e["left"], depth=e["depth"]))
print("%d phases with histories from %s" % (len(ph), " ".join(sys.argv[2:])))
print("| rule | phase it | children [e] | total | x recorded | lost fathoms | stopped early | censored (would run on) |")
print("|---|---|---|---|---|---|---|---|")
ref = None
for name, rule in RULES:
    pit = kid = lost = early = cens = 0
    for p in ph:
        h = p["hist"]; stop = None
        for k in range(len(h)):
            if h[k][1] <= 0: break
            if rule(p, h, k): stop = k; break
        if stop is not None and h[stop][0] < p["its"]:
            pit += h[stop][0]; kid += p["m"] * child(h[stop][1] / UNIT); early += 1; lost += p["fath"]
        else:
            pit += p["its"]
            if not p["fath"]: kid += p["m"] * child(p["left"] / UNIT); cens += stop is None
    tot = pit + kid
    if ref is None: ref = tot
    print("| %s | %d | %d | %d | %.2f | %d | %d | %d |" % (name, pit, kid, tot, ref / tot, lost, early, cens))
