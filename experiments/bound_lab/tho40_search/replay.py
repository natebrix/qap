"""What the measured cut phase would do to a production (no-cut) node trace: a replay model [e].

  python3 replay.py [--trace T.gz] [--zone Z ...] [--only SPEC] [--sf S] [--carry all|none]

The production trace (infra/phase1/runs/tho40ab-32/nodes.txt.gz: the A/B subtree -F 26:10,8:0, tho40.par's
settings, no cut phase) has, for every branched node, the shortfall s at its warm give-up and the CPU of
everything under it. The model walks it top-down in production CPU seconds:

  node fathomed in production (D)       its own CPU
  branched, s > zone                    its own CPU + its children (each by this table)
  branched, s <= zone, s <= SF          its own CPU + a phase that fathoms: ITS_F(s) iterations
  branched, s <= zone, s > SF           its own CPU + a phase that misses: ITS_M iterations, then its children:
       --carry all   every child fathoms at its pooled first solve (POOLED x the mean CPU of a fathoming
                     node at that depth): what the two missed root phases here did (t3, t4 z030: 73 of 73)
       --carry none  every child goes on as if its parent had run no phase (by this table)

SF, ITS_F, ITS_M and POOLED are this directory's measurements (phases.py; README.md, "The phases").
A phase iteration is charged S_PER_IT production CPU seconds at m = 36, scaled by (m/36)^3: 0.118 s here
(t2 z030: 2,070 iterations, 244 s) over 0.70, this box's CPU per matched node against production's (calib.py).
"""
import sys, argparse, collections
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
import ntrace

INC = 240516.0
S_PER_IT = 0.118 / 0.70
# phase iterations to fathom from shortfall s [m]: warm children 88:44 185:59 232:73 236:78 411:140 823:276
# 1247:526 1443:773; cold roots 385:130 957:220 3193:2070
PTS = [(0, 30), (250, 80), (400, 140), (800, 270), (1250, 530), (1450, 770), (3200, 2070), (5000, 3400)]
ITS_M = 2600      # a miss: 1,810 (stalled, t4) and 3,400 (the limit, t3)
POOLED = 1.5      # a pooled child's first solve against an unpooled fathoming node's (t3, t4 z030: 19-45 it against 15-25)

def its_f(s):
    for (a, b), (c, d) in zip(PTS, PTS[1:]):
        if s <= c: return b + (d - b) * (s - a) / (c - a)
    return PTS[-1][1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", default="/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
    ap.add_argument("--zone", type=float, nargs="*", default=[0.003, 0.01, 0.03])
    ap.add_argument("--only", help="the subtree of this -F spec only")
    ap.add_argument("--sf", type=float, nargs="*", default=[3000.0], help="a phase fathoms from at most this shortfall")
    ap.add_argument("--carry", nargs="*", default=["all", "none"])
    ap.add_argument("--by-depth", action="store_true")
    a = ap.parse_args()
    t = ntrace.load(a.trace)
    roots = t.roots
    if a.only:
        fx = tuple(a.only.split(","))
        roots = [i for i in range(len(t)) if t.fixes[i] == fx]
    order = sorted(range(len(t)), key=lambda i: -t.depth[i])
    dn = collections.Counter(); dc = collections.Counter()
    for i in range(len(t)):
        if t.outcome[i] != "B": dn[t.depth[i]] += 1; dc[t.depth[i]] += t.cpu[i]
    dmean = {d: dc[d] / dn[d] for d in dn}
    base_n = sum(t.size[r] for r in roots); base_c = sum(t.sub_cpu[r] for r in roots)
    print("production: %d nodes, %.0f CPU s" % (base_n, base_c))
    print("zone   short   SF  carry |  nodes  (x fewer) | CPU s  (x less) | phases: fathoming, missing | share of CPU in phases")
    for Z in a.zone:
        for SF in a.sf:
            for carry in a.carry:
                En = [0.0] * len(t); Ec = [0.0] * len(t); Pf = [0.0] * len(t); Pm = [0.0] * len(t); Pc = [0.0] * len(t)
                for i in order:
                    n = 1.0; c = t.cpu[i]; pf = pm = pc = 0.0
                    if t.outcome[i] == "B":
                        s = t.cutoff[i] - t.dnn[i]; m = 40 - t.depth[i]; ch = t.children[i]
                        unit = S_PER_IT * (m / 36.0) ** 3
                        if s > Z * INC:
                            n += sum(En[k] for k in ch); c += sum(Ec[k] for k in ch)
                            pf = sum(Pf[k] for k in ch); pm = sum(Pm[k] for k in ch); pc = sum(Pc[k] for k in ch)
                        elif s <= SF:
                            c += unit * its_f(s); pf = 1; pc = unit * its_f(s)
                        else:
                            c += unit * ITS_M; pm = 1; pc = unit * ITS_M
                            if carry == "all":
                                n += len(ch); c += len(ch) * POOLED * dmean.get(t.depth[i] + 1, 2.0)
                            else:
                                n += sum(En[k] for k in ch); c += sum(Ec[k] for k in ch)
                                pf += sum(Pf[k] for k in ch); pm += sum(Pm[k] for k in ch); pc += sum(Pc[k] for k in ch)
                    En[i] = n; Ec[i] = c; Pf[i] = pf; Pm[i] = pm; Pc[i] = pc
                mn = sum(En[r] for r in roots); mc = sum(Ec[r] for r in roots)
                print("%.3f %6.0f %5.0f %5s | %6.0f (x%5.1f) | %6.0f (x%4.2f) | %4.0f, %3.0f | %.0f%%" % (
                    Z, Z * INC, SF, carry, mn, base_n / mn, mc, base_c / mc, sum(Pf[r] for r in roots), sum(Pm[r] for r in roots),
                    100 * sum(Pc[r] for r in roots) / mc))

if __name__ == "__main__":
    main()
