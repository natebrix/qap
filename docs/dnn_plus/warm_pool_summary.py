#!/usr/bin/env python3
"""Tables for warm_pool.sh: python3 docs/dnn_plus/warm_pool_summary.py [OUT]
For each case: the children every arm of an experiment bounded, how many each arm fathomed,
decisions that differ from the experiment's first arm, and iterations and CPU seconds
(a) on the children fathomed in every arm and (b) on all children."""
import os, re, sys
H = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{H}/runs_warm"
EXPS = [("A: warm starts", ["cold", "z", "zy"]), ("B: inherited cuts", ["nopool", "pool"])]

def read(path):
    """fixes -> dict(depth, lb, fath, iters, cpu, warm, pool), one per node, and the parent's key"""
    nodes = {}
    for l in open(path):
        m = re.search(r"depth (\d+) fixes ?(.*?) ?\| first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+) warm (\d+) pool (\d+)", l)
        if not m: continue
        cut = re.search(r"\| cut (.*?) \| cutoff", l).group(1)
        f = int(re.search(r"cutoff \S+ fathomed (\d)", l).group(1))
        nodes[m.group(2)] = dict(depth=int(m.group(1)), lb=float(m.group(3)), status=m.group(4), iters=int(m.group(5)),
                                 cpu=float(m.group(6)), warm=int(m.group(7)), pool=int(m.group(8)), fath=f, cut=cut)
    return nodes

for case in sorted(os.listdir(OUT)):
    d = f"{OUT}/{case}"
    if not os.path.isdir(d): continue
    print(f"== {case}  ({open(d + '/info.txt').read().strip() if os.path.exists(d + '/info.txt') else ''})")
    for name, arms in EXPS:
        t = {a: read(f"{d}/{a}.trace") for a in arms if os.path.exists(f"{d}/{a}.trace")}
        if len(t) < 2: continue
        arms = [a for a in arms if a in t]
        dmin = min(n["depth"] for n in t[arms[0]].values())
        kids = {a: {k: n for k, n in t[a].items() if n["depth"] > dmin} for a in arms}
        par = {a: [n for n in t[a].values() if n["depth"] == dmin] for a in arms}
        common = sorted(set.intersection(*[set(kids[a]) for a in arms]))
        allf = [k for k in common if all(kids[a][k]["fath"] for a in arms)]
        print(f"  {name}: children {', '.join(f'{a} {len(kids[a])}' for a in arms)}; in every arm {len(common)}; fathomed in every arm {len(allf)}")
        for a in arms:
            p = par[a][0] if par[a] else None
            if p: print(f"    parent {a:7s} first lb {p['lb']:.2f} iters {p['iters']} fathomed {p['fath']} cut [{p['cut'][:70]}]")
        print(f"    {'arm':7s} {'fathomed':>8s} {'changed':>7s} {'it(shared fathoms)':>18s} {'cpu':>8s} {'it(all)':>9s} {'cpu':>8s} {'warm':>5s} {'pooled':>6s} {'mean lb, never fathomed':>24s}")
        for a in arms:
            K = kids[a]
            nf = [k for k in common if not any(kids[b][k]["fath"] for b in arms)]
            ch = sum(K[k]["fath"] != kids[arms[0]][k]["fath"] for k in common)
            print(f"    {a:7s} {sum(K[k]['fath'] for k in common):8d} {ch:7d} {sum(K[k]['iters'] for k in allf):18d} {sum(K[k]['cpu'] for k in allf):8.0f}"
                  f" {sum(K[k]['iters'] for k in common):9d} {sum(K[k]['cpu'] for k in common):8.0f}"
                  f" {sum(K[k]['warm'] > 0 for k in common):5d} {sum(K[k]['pool'] > 0 for k in common):6d}"
                  f" {(sum(K[k]['lb'] for k in nf) / len(nf) if nf else 0):18.2f} ({len(nf)})")
