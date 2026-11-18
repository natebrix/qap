"""Replay Knuth-estimator designs on a finished branch-and-bound tree (docs/estimates_plan.md, step 1).

A finished run's node trace has every node: its fixes, its bounds (QPB, the DNN first solve, the cut
phase), how it ended and its CPU. The replay takes the estimator's dives on that tree instead of
solving nodes, so each design's estimates can be scored against the known totals, at the CPU a real
run of the design would spend (the recorded CPU of every node it would bound).

    python3 -I replay.py TRACE [--k 2,3] [--budget 0.01,0.03,0.1] [--reps 200] [--designs ...]

TRACE: a node trace (QAP_NODE_TRACE, `qapctl export`'s nodes.txt[.gz]) or a directory of the website's
tree shards (tree/<run>/data/*.json).

Designs (all start from the exact top: every node at depth < k counted as it is, every depth-k node
bounded once; a depth-k node that does not branch is exact):
  y2000   the thesis's and today's estimate.C: each dive starts at a depth-k node drawn by g^q, and takes
          one child, drawn by g^q (g the relative gap of the child's own QPB bound against the root's, as
          QAP_EST_CHILD_QPB=1; 0.1 uniform mixed in, est.xbar_mix's default)
  strat   the same dives, but every branched depth-k node is its own stratum with its own dives
  surv    stratified, and at every node of a dive all children are bounded (the search bounds them all
          anyway); the fathomed ones are counted exactly and the dive goes on into a child that
          branched, drawn uniformly (surv) or by its shortfall^q (survq: cutoff - its best bound; 0.1
          uniform mixed in). Unbiased (Purdom 1978, partial backtracking)
Allocation within the stratified designs: a pilot of 2 dives per stratum, then by s_h / sqrt(c_h)
(Neyman at fixed cost), s_h floored at the stratum's mean.
Cost: a node bounded by several dives of one run is paid once (the memo of docs/estimates_plan.md).
"""
import argparse, bisect, glob, gzip, json, math, os, random, statistics, sys
from collections import defaultdict


def load_trace(path):
    op = gzip.open if path.endswith(".gz") else open
    nodes = {}
    with op(path, "rt") as f:
        for line in f:
            t = line.split()
            if len(t) < 12 or t[0] != "node": continue
            fx = "" if t[2] == "-" else t[2]
            if fx in nodes: continue                      # a unit's retry: keep the first
            v = lambda s: float("nan") if s == "-" else float(s)
            nodes[fx] = dict(depth=int(t[1]), qpb=v(t[3]), dnn=v(t[4]), cut=v(t[5]), cutoff=v(t[6]),
                             out=t[7], cpu=float(t[11]))
    return nodes


def load_shards(d):
    rows = {}
    for p in glob.glob(os.path.join(d, "*.json")):
        if p.endswith("meta.json"): continue
        s = json.load(open(p)); c = s["columns"]
        for r in s["rows"]:
            x = dict(zip(c, r)); rows[x["id"]] = x
    nodes, key = {}, {}
    for i in sorted(rows, key=lambda i: rows[i]["depth"]):
        x = rows[i]
        k = "" if x["parent"] in (-1, None) else (key[x["parent"]] + "," if key[x["parent"]] else "") + x["fix"]
        key[i] = k
        f = lambda v: float("nan") if v is None else float(v)
        nodes[k] = dict(depth=x["depth"], qpb=f(x["qpb"]), dnn=f(x["dnn"]), cut=f(x["cut"]),
                        cutoff=f(x["cutoff"]), out=x["outcome"], cpu=f(x["cpu"]) if x["cpu"] is not None else 0.0)
    return nodes


class Tree:
    def __init__(self, nodes):
        keys = sorted(nodes, key=lambda k: (nodes[k]["depth"], k))
        self.id = {k: i for i, k in enumerate(keys)}
        n = len(keys)
        self.depth = [nodes[k]["depth"] for k in keys]
        self.out = [nodes[k]["out"] for k in keys]
        self.cpu = [nodes[k]["cpu"] for k in keys]
        self.cnt = [0 if nodes[k]["out"] == "P" else 1 for k in keys]
        self.qpb = [nodes[k]["qpb"] for k in keys]
        best = []
        for k in keys:
            x = nodes[k]; b = x["cut"] if not math.isnan(x["cut"]) else (x["dnn"] if not math.isnan(x["dnn"]) else x["qpb"])
            best.append(b)
        self.best = best
        self.cutoff = [nodes[k]["cutoff"] for k in keys]
        self.kids = [[] for _ in range(n)]
        orphans = 0
        for k in keys:
            if not k: continue
            p = k.rsplit(",", 1)[0] if "," in k else ""
            if p in self.id: self.kids[self.id[p]].append(self.id[k])
            else: orphans += 1
        self.orphans = orphans
        self.root = self.id[""]
        self.opt = max(c for c in self.cutoff if not math.isnan(c)) + 0.95   # the first cutoff: incumbent - 1 + 0.05
        self.root_qpb = self.qpb[self.root]
        # subtree totals
        self.N = [0] * n; self.C = [0.0] * n; self.maxd = [0] * n
        for i in reversed(range(n)):
            self.N[i] = self.cnt[i] + sum(self.N[j] for j in self.kids[i])
            self.C[i] = self.cpu[i] + sum(self.C[j] for j in self.kids[i])
            self.maxd[i] = max([self.depth[i]] + [self.maxd[j] for j in self.kids[i]])
        self.branched = [self.out[i] == "B" and bool(self.kids[i]) for i in range(n)]
        self._prep = {}

    def g(self, b):
        if not (b < self.opt): return 0.0
        return (self.opt - b) / (self.opt - self.root_qpb)

    def knuth_probs(self, i, q, mix=0.1):
        key = ("k", i, q)
        if key not in self._prep:
            ks = self.kids[i]
            w = [self.g(self.qpb[j]) ** q if self.g(self.qpb[j]) > 0 else 0.0 for j in ks]
            s = sum(w)
            p = [((1 - mix) * x / s + mix / len(ks)) if s > 0 else 1 / len(ks) for x in w]
            self._prep[key] = (ks, p, cum(p))
        return self._prep[key]

    def surv_probs(self, i, q, mix=0.1):
        key = ("s", i, q)
        if key not in self._prep:
            ks = self.kids[i]
            surv = [j for j in ks if self.branched[j]]
            fn = sum(self.cnt[j] for j in ks if not self.branched[j])
            fc = sum(self.cpu[j] for j in ks if not self.branched[j])
            if surv:
                if q == 0: p = [1 / len(surv)] * len(surv)
                else:
                    w = [max(self.cutoff[j] - self.best[j], 0.0) ** q for j in surv]; s = sum(w)
                    p = [((1 - mix) * x / s + mix / len(surv)) if s > 0 else 1 / len(surv) for x in w]
            else: p = []
            self._prep[key] = (surv, p, cum(p), fn, fc)
        return self._prep[key]


def cum(p):
    c, s = [], 0.0
    for x in p: s += x; c.append(s)
    return c


def draw(rng, cp):
    u = rng.random() * cp[-1]
    return min(bisect.bisect_left(cp, u), len(cp) - 1)


def dive(T, start, d0, mode, q, rng, paid):
    """one dive from node start with weight d0; returns (nodes, cpu) estimate and the CPU newly paid"""
    en = ec = 0.0; cost = 0.0; i = start; d = d0; last = T.depth[start]
    while True:
        if i not in paid: paid.add(i); cost += T.cpu[i]
        en += d * T.cnt[i]; ec += d * T.cpu[i]; last = T.depth[i]
        if not T.branched[i]: break
        if mode == "knuth":
            ks, p, cp = T.knuth_probs(i, q)
            k = draw(rng, cp); i, d = ks[k], d / p[k]
        else:
            for j in T.kids[i]:
                if j not in paid: paid.add(j); cost += T.cpu[j]
            surv, p, cp, fn, fc = T.surv_probs(i, q)
            en += d * fn; ec += d * fc
            if not surv:
                if T.kids[i]: last = T.depth[i] + 1
                break
            k = draw(rng, cp); j = surv[k]; d = d / p[k]
            # j's own count and cpu are added at the top of the loop, with its weight d
            i = j
            # j is already paid (bounded as a child)
    return en, ec, cost, last


MAXDIVES = 20_000   # dives the memo makes free (every node on their paths paid) end a replication here


def run_design(T, design, k, budget_cpu, rng, q):
    """one replication: returns (est nodes, est cpu, cpu spent, deepest depth reached)"""
    n = len(T.depth)
    top = [i for i in range(n) if T.depth[i] < k]
    strata = [i for i in range(n) if T.depth[i] == k]
    paid = set(top) | set(strata)
    spent = sum(T.cpu[i] for i in paid)
    base_n = sum(T.cnt[i] for i in top); base_c = sum(T.cpu[i] for i in top)
    exact = [s for s in strata if not T.branched[s]]
    base_n += sum(T.cnt[s] for s in exact); base_c += sum(T.cpu[s] for s in exact)
    br = [s for s in strata if T.branched[s]]
    if not br: return base_n, base_c, spent, k
    mode = "knuth" if design in ("y2000", "strat") else "surv"
    deep = k
    def one(s, d0):
        nonlocal spent, deep
        en, ec, c, last = dive(T, s, d0, mode, q, rng, paid)
        spent += c; deep = max(deep, last)
        return en, ec, c
    if design == "y2000":
        w = [T.g(T.qpb[s]) ** q if T.g(T.qpb[s]) > 0 else 0.0 for s in br]; sw = sum(w)
        p = [(0.9 * x / sw + 0.1 / len(br)) if sw > 0 else 1 / len(br) for x in w]; cp = cum(p)
        sn = sc = 0.0; m = 0
        while (spent < budget_cpu or m < 2) and m < MAXDIVES:
            j = draw(rng, cp); s = br[j]
            en, ec, _ = one(s, 1.0 / p[j]); sn += en; sc += ec; m += 1
        return base_n + sn / m, base_c + sc / m, spent, deep
    # stratified: pilot of 2, then Neyman at fixed cost (running sums; a stratum's cost per dive is the
    # mean of its last dives', floored at 1 CPU s, since the memo drives it toward 0)
    st = {s: [0, 0.0, 0.0, 0.0, 0.0] for s in br}      # dives, sum nodes, sum cpu, sum cpu^2, cost EMA
    def add(s):
        en, ec, c = one(s, 1.0); x = st[s]
        x[0] += 1; x[1] += en; x[2] += ec; x[3] += ec * ec
        x[4] = c if x[0] == 1 else 0.7 * x[4] + 0.3 * c
    for s in br:
        add(s); add(s)
    total = 2 * len(br)
    while spent < budget_cpu and total < MAXDIVES:
        wts = []
        for s in br:
            n_, s1, s2 = st[s][0], st[s][2], st[s][3]
            mu = s1 / n_; sd = math.sqrt(max(s2 / n_ - mu * mu, 0.0)); sd = max(sd, mu)
            wts.append(sd / math.sqrt(max(st[s][4], 1.0)) / n_)
        if sum(wts) <= 0: break
        cpw = cum(wts)
        for _ in range(max(1, len(br) // 4)):
            add(br[draw(rng, cpw)]); total += 1
            if spent >= budget_cpu: break
    en = base_n + sum(st[s][1] / st[s][0] for s in br)
    ec = base_c + sum(st[s][2] / st[s][0] for s in br)
    return en, ec, spent, deep


def run_split(T, k, budget_cpu, rng, q, fresh, structure=0.5, n0=2, by="var"):
    """adaptive stratification with survivor dives: from the exact top at depth k, split (bound all children
    of, as one more breadth-first level) the stratum with the largest variance contribution s_h^2 / n_h,
    until `structure` of the budget is spent; then Neyman dives to the budget. fresh: the final estimate
    uses only the dives taken after the structure was fixed (unbiased), else every dive of a stratum."""
    n = len(T.depth)
    top = [i for i in range(n) if T.depth[i] < k]
    strata = [i for i in range(n) if T.depth[i] == k]
    paid = set(top) | set(strata)
    spent = sum(T.cpu[i] for i in paid)
    base_n = sum(T.cnt[i] for i in top) + sum(T.cnt[s] for s in strata if not T.branched[s])
    base_c = sum(T.cpu[i] for i in top) + sum(T.cpu[s] for s in strata if not T.branched[s])
    H = [s for s in strata if T.branched[s]]
    st = {}
    deep = k
    def add(s):
        nonlocal spent, deep
        en, ec, c, last = dive(T, s, 1.0, "surv", q, rng, paid)
        spent += c; deep = max(deep, last)
        x = st.setdefault(s, [0, 0.0, 0.0, 0.0, 0.0])
        x[0] += 1; x[1] += en; x[2] += ec; x[3] += ec * ec
        x[4] = c if x[0] == 1 else 0.7 * x[4] + 0.3 * c
    def var(s):
        x = st[s]; mu = x[2] / x[0]; return max(x[3] / x[0] - mu * mu, 0.0) / x[0]
    for s in H:
        for _ in range(n0): add(s)
    key = {"var": var, "mean": lambda s: st[s][2] / st[s][0], "short": lambda s: T.cutoff[s] - T.best[s],
           "bfs": lambda s: 1000 - T.depth[s] + 1e-9 * (T.cutoff[s] - T.best[s])}[by]   # bfs: the shallowest first, as a solve's frontier grows
    while spent < structure * budget_cpu and H:
        h = max(H, key=key)
        if key(h) <= 0: break
        H.remove(h); st.pop(h, None)
        base_n += T.cnt[h]; base_c += T.cpu[h]
        for j in T.kids[h]:
            if j not in paid: paid.add(j); spent += T.cpu[j]
            if T.branched[j]:
                H.append(j)
                for _ in range(n0): add(j)
            else:
                base_n += T.cnt[j]; base_c += T.cpu[j]
    if fresh:
        st.clear()
        for s in H:
            for _ in range(n0): add(s)
    total = sum(st[s][0] for s in H)
    while spent < budget_cpu and total < MAXDIVES and H:
        wts = []
        for s in H:
            x = st[s]; mu = x[2] / x[0]; sd = math.sqrt(max(x[3] / x[0] - mu * mu, 0.0)); sd = max(sd, mu)
            wts.append(sd / math.sqrt(max(x[4], 1.0)) / x[0])
        if sum(wts) <= 0: break
        cpw = cum(wts)
        for _ in range(max(1, len(H) // 4)):
            add(H[draw(rng, cpw)]); total += 1
            if spent >= budget_cpu: break
    en = base_n + sum(st[s][1] / st[s][0] for s in H)
    ec = base_c + sum(st[s][2] / st[s][0] for s in H)
    return en, ec, spent, deep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trace")
    ap.add_argument("--k", default="2,3")
    ap.add_argument("--budget", default="0.01,0.03,0.1", help="the estimate's CPU as a fraction of the tree's")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--designs", default="y2000:1.5,strat:1.5,surv:0,survq:1")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    nodes = load_shards(a.trace) if os.path.isdir(a.trace) else load_trace(a.trace)
    T = Tree(nodes)
    N, C = T.N[T.root], T.C[T.root]
    print(f"# {a.trace}: {N:,} nodes, {C/3600:,.1f} CPU h, deepest {T.maxd[T.root]}, {T.orphans} orphans")
    print("| k | design | q | budget | median est/true nodes | 10-90% | median est/true CPU | 10-90% | within 20% (CPU) | within 2x (CPU) | CPU spent / true | deepest reached (median) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    rng = random.Random(a.seed)
    for k in map(int, a.k.split(",")):
        for dq in a.designs.split(","):
            design, q = dq.split(":"); q = float(q)
            for b in map(float, a.budget.split(",")):
                rn, rc, sp, dp = [], [], [], []
                for _ in range(a.reps):
                    if design.startswith("split"):
                        # split[fresh][-var|-mean|-short][-sN: the structure's share of the budget, %]
                        parts = design.split("-"); by = next((x for x in parts[1:] if x in ("var", "mean", "short", "bfs")), "var")
                        sh = next((int(x[1:]) / 100 for x in parts[1:] if x.startswith("s") and x[1:].isdigit()), 0.5)
                        en, ec, spent, dd = run_split(T, k, b * C, rng, q, parts[0] == "splitfresh", structure=sh, by=by)
                    else:
                        en, ec, spent, dd = run_design(T, "surv" if design == "survq" else design, k, b * C, rng, q)
                    rn.append(en / N); rc.append(ec / C); sp.append(spent / C); dp.append(dd)
                rn.sort(); rc.sort()
                pct = lambda xs, f: xs[min(int(f * len(xs)), len(xs) - 1)]
                w20 = sum(0.8 <= x <= 1.2 for x in rc) / len(rc); w2 = sum(0.5 <= x <= 2 for x in rc) / len(rc)
                print(f"| {k} | {design} | {q:g} | {b:g} | {pct(rn,.5):.2f} | {pct(rn,.1):.2f}-{pct(rn,.9):.2f} | "
                      f"{pct(rc,.5):.2f} | {pct(rc,.1):.2f}-{pct(rc,.9):.2f} | {w20:.0%} | {w2:.0%} | {statistics.fmean(sp):.3f} | {statistics.median(dp):.0f} |", flush=True)


if __name__ == "__main__":
    main()
