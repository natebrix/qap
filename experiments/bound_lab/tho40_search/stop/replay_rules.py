"""Replay candidate stop rules on the recorded phase histories (hist.py), offline.

A rule sees what the solver's stop() sees at each check (every 10 iterations): the iteration and the
running best. It returns the iteration at which it would stop the phase, or None. A phase the rule does
not stop before its recorded end keeps its recorded outcome (so a rule more patient than the recorded
one, cut.stall 3 x cut.iters 200, is scored as the record: censored, counted in column 'cens').

Cost of a phase, in iterations:  t_stop  (+ C(gap at t_stop) if it did not fathom)
C(g) is what the node's children then cost, HANDOVER below: m x (N0 + A g + B g^2), from the misses whose
children ran complete (printed first, model against measured). --kappa K weights a phase iteration K x a
first-solve iteration (CPU: 1.25); default 1 = iterations, the metric of the README.

  python3 replay_rules.py [--kappa 1.25] [--runs w18_zall,w38_z030,...] [--detail RULE]
"""
import sys, argparse, math
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from hist import load_phases

# children's cost after a hand-over at gap g, per child: 20 iterations when the parent ended at the cutoff
# (t3 z030: 19.1; w18's 2:7: 20.6), rising with the gap as more of them need a phase of their own.
# g in 0.1% of the cutoff (UNIT = 240.514 units), as the solver key cut.handover = N0:A:B takes it.
UNIT = 240.514
N0, A, B = 20.0, 2.0, 3.0
def handover(g, m, n0=N0, a=A, b=B):
    u = g / UNIT
    return m * (n0 + a * u + b * u * u)

def gain_over(h, k, W):
    """running-best gain over the last W iterations at check k (h: [(it, gap)]); gap0 before the first check"""
    it, g = h[k]
    g0 = h[0][1] if False else None
    prev = None
    for jt, gj in h:
        if jt > it - W: break
        prev = gj
    return prev

def rule_stall(W, f=3.0, tmin=None):
    def r(p, h, k):
        it, g = h[k]
        if it < (tmin or W): return False
        g0 = p["short0"]
        for jt, gj in h:
            if jt > it - W: break
            g0 = gj
        return g > f * max(0.0, g0 - g)
    return r

def rule_cap(T, base):
    def r(p, h, k): return h[k][0] >= T or base(p, h, k)
    return r

def rule_giveup(pairs, base):
    def r(p, h, k):
        it, g = h[k]
        return any(it >= T and g > G for T, G in pairs) or base(p, h, k)
    return r

def rule_handover(W=200, phi=1.5, n0=N0, a=A, b=B):
    """cut.handover = n0:a:b:phi:W. Stop when the iterations the phase still needs at its present rate,
    phi x gap / rate, exceed what handing over to the children would cost, m x (n0 + a g + b g^2);
    rate = gain of the last W iterations / W"""
    def r(p, h, k):
        it, g = h[k]
        if it < W: return False
        g0 = p["short0"]
        for jt, gj in h:
            if jt > it - W: break
            g0 = gj
        rate = max(0.0, g0 - g) / W
        return phi * g > rate * handover(g, p["m"], n0, a, b)
    return r

base = rule_stall(600)
RULES = [
    ("recorded: stall 3 x iters 200", base),
    ("stall 2 x iters 200 (W 400)", rule_stall(400)),
    ("stall 3 x iters 100 (W 300)", rule_stall(300)),
    ("stall 2 x iters 100 (W 200)", rule_stall(200)),
    ("cut.cp_iters 1500", rule_cap(1500, base)),
    ("cut.cp_iters 1000", rule_cap(1000, base)),
    ("cut.cp_iters 600", rule_cap(600, base)),
    ("cut.giveup 600:200", rule_giveup([(600, 200)], base)),
    ("cut.giveup 400:400,800:150", rule_giveup([(400, 400), (800, 150)], base)),
    ("cut.giveup 300:250,600:100", rule_giveup([(300, 250), (600, 100)], base)),
    ("rate: gap/rate > 20 m (W 200)", rule_handover(200, 1.0, 20, 0, 0)),
    ("rate: gap/rate > 40 m (W 200)", rule_handover(200, 1.0, 40, 0, 0)),
    ("handover 20:2:3 phi 1.5 W 200", rule_handover(200, 1.5)),
    ("handover 20:2:3 phi 1 W 200", rule_handover(200, 1.0)),
    ("handover 20:2:3 phi 2 W 200", rule_handover(200, 2.0)),
    ("handover 20:2:3 phi 1.5 W 100", rule_handover(100, 1.5)),
    ("handover 20:2:3 phi 1.5 W 300", rule_handover(300, 1.5)),
    ("handover 30:2:3 phi 1.5 W 200", rule_handover(200, 1.5, 30)),
    ("handover 20:2:1.5 phi 1.5 W 200", rule_handover(200, 1.5, 20, 2, 1.5)),
    ("handover 20:2:6 phi 1.5 W 200", rule_handover(200, 1.5, 20, 2, 6)),
]

def play(p, rule):
    h = p["hist"]
    for k in range(len(h)):
        it, g = h[k]
        if g < 0: break                      # fathomed here
        if it >= p["its"] and p["status"] != "fathomed": break
        if rule(p, h, k): return it, g, False
    return p["its"], p["left"], True         # ran to its recorded end

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--runs")
    ap.add_argument("--detail")
    a = ap.parse_args()
    P = [p for p in load_phases(a.runs.split(",") if a.runs else None) if not p["path"]]
    if not a.runs: P = [p for p in P if "zallh" not in p["run"]]   # the phases recorded under the stall rule only
    print("the hand-over model against the misses whose children ran: gap left, m, children measured, model")
    for p in P:
        if p["after"] is not None and p["status"] != "fathomed" and p["after"] > 0:
            print("   %-10s d%d left %7.1f  m %d  measured %6d  model %6.0f" % (p["run"], p["depth"], p["left"], p["m"], p["after"], handover(p["left"], p["m"])))
    near = [p for p in P if p["short0"] <= 2000]; far = [p for p in P if p["short0"] > 2000]
    print("\n%d phases (%s): %d that start within 2,000 of the cutoff (%d fathomed as recorded), %d farther (%d fathomed)" % (
        len(P), a.runs or "all runs with a history", len(near), sum(p["status"] == "fathomed" for p in near), len(far), sum(p["status"] == "fathomed" for p in far)))
    print("cost = kappa %.2f x phase iterations + the modelled hand-over after a stop. 'lost': recorded fathoms the rule stops first;" % a.kappa)
    print("'more': recorded stalls the rule would have let run on (scored as recorded: their extension is not in the data)")
    print("%-32s | near: phase it  handover  total  lost | far: phase it  handover  total  lost more | all" % "rule")
    ref = None
    for name, rule in RULES:
        res = []
        for S in (near, far):
            tp = th = 0.0; lost = more = 0
            for p in S:
                t, g, end = play(p, rule)
                c = 0.0
                if end:
                    if p["status"] != "fathomed":
                        c = handover(p["left"], p["m"])
                        if p["status"] == "stalled" and not rule(p, p["hist"], len(p["hist"]) - 1): more += 1
                else:
                    c = handover(g, p["m"]); lost += p["status"] == "fathomed"
                tp += t; th += c
                if a.detail and a.detail in name:
                    print("      %-10s d%d %s start %7.1f: stop at %4d gap %7.1f (%s; recorded %s %d it) cost %6.0f against recorded %6.0f" % (
                        p["run"], p["depth"], "pool" if p["pooled"] else "    ", p["short0"], t, g, "recorded end" if end else "rule", p["status"], p["its"],
                        a.kappa * t + c, a.kappa * p["its"] + (handover(p["left"], p["m"]) if p["status"] != "fathomed" else 0)))
            res.append((tp, th, a.kappa * tp + th, lost, more))
        tot = res[0][2] + res[1][2]
        if ref is None: ref = (res[0][2], res[1][2], tot)
        print("%-32s | %8.0f %8.0f %6.0f (x%.2f) %2d | %8.0f %8.0f %6.0f (x%.2f) %2d %2d | %6.0f (x%.2f)" % (
            name, res[0][0], res[0][1], res[0][2], ref[0] / res[0][2], res[0][3], res[1][0], res[1][1], res[1][2], ref[1] / res[1][2], res[1][3], res[1][4], tot, ref[2] / tot))

if __name__ == "__main__":
    main()
