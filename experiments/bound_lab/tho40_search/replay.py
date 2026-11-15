"""What the measured cut phase would do to a production (no-cut) node trace: a replay model [e].

  python3 replay.py [--trace T.gz] [--zone Z ...] [--only SPEC]

The production trace has, for every branched node, the shortfall at its (warm) give-up and the CPU of
everything under it. The model walks it top-down. A branched node whose shortfall is within the zone
(Z x 240,516) runs a phase:
  - it fathoms with probability p(s), s the shortfall, and costs it_f(s) phase iterations;
  - else it costs it_m(s) iterations and its children go on as in production (the pool the children
    would inherit is ignored: conservative for the cut arm).
p, it_f and it_m are the measured phases of this directory's runs (phases.py), binned by shortfall
(MODEL below; edit when more phases are measured). An iteration of a phase is charged
PHASE_S_PER_IT_PROD production CPU seconds at m = 36, scaled by (m/36)^3.
The expected cost is E[node] = cpu(node) + [in zone] (phase cost + (1 - p) sum E[children]),
and nodes E_n[node] = 1 + (1 - p) sum E_n[children].
"""
import sys, os, argparse
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
import ntrace

INC = 240516.0
# shortfall upper edge, p(fathom), iterations when it fathoms, iterations when it misses  [m: phases.py]
MODEL = [
    (500, 1.0, 130, 600),
    (1000, 1.0, 220, 600),
    (1500, 1.0, 526, 800),
    (1e9, 0.0, 0, 800),
]
PHASE_S_PER_IT_PROD = 0.17   # production CPU s per phase iteration at m = 36 (calib.py: 0.11 s here / 0.65)

def model(s):
    for hi, p, itf, itm in MODEL:
        if s <= hi: return p, itf, itm
    return 0.0, 0, 800

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", default="/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
    ap.add_argument("--zone", type=float, nargs="*", default=[0.003, 0.01, 0.03])
    ap.add_argument("--only", help="the subtree of this -F spec only")
    ap.add_argument("--root-phase", type=int, default=1, help="0: no phase at the subtree's root (a cold root outside the zone)")
    a = ap.parse_args()
    t = ntrace.load(a.trace)
    roots = t.roots
    if a.only:
        fx = tuple(a.only.split(","))
        roots = [i for i in range(len(t)) if t.fixes[i] == fx]
    order = sorted(range(len(t)), key=lambda i: -t.depth[i])
    base_n = sum(t.size[r] for r in roots); base_c = sum(t.sub_cpu[r] for r in roots)
    print("production: %d nodes, %.0f CPU s" % (base_n, base_c))
    for Z in a.zone:
        En = [0.0] * len(t); Ec = [0.0] * len(t); phases = 0.0; pc = 0.0
        reach = [0.0] * len(t)
        for i in order:
            n = 1.0; c = t.cpu[i]
            if t.outcome[i] == "B":
                s = t.cutoff[i] - t.dnn[i]
                kn = sum(En[k] for k in t.children[i]); kc = sum(Ec[k] for k in t.children[i])
                if s <= Z * INC and not (i in roots and not a.root_phase):
                    p, itf, itm = model(s)
                    m = 40 - t.depth[i]
                    cost = PHASE_S_PER_IT_PROD * (m / 36.0) ** 3 * (p * itf + (1 - p) * itm)
                    n += (1 - p) * kn; c += cost + (1 - p) * kc
                else:
                    n += kn; c += kc
            En[i] = n; Ec[i] = c
        mn = sum(En[r] for r in roots); mc = sum(Ec[r] for r in roots)
        print("zone %.4f (%5.0f short): nodes %8.0f (x%.2f fewer), CPU %8.0f (x%.2f less)" % (Z, Z * INC, mn, base_n / mn, mc, base_c / mc))

if __name__ == "__main__":
    main()
