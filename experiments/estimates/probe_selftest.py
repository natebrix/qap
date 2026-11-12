#!/usr/bin/env python3
"""probe.py's estimator on synthetic trees (no solver): the mean of many probes
against the true size, one probe's spread, and how often +-1.96 SE holds the truth.

  probe_selftest.py [--beta 1.8] [--k 12] [--reps 60] [--s0 2.6]

A node s > 0 short has 3-12 children, each 0.4-1.6 shorter (a node at s <= 0 is
fathomed) and costs 10 + 50 s u CPU s (u uniform), all fixed by the node's id, so
the tree is the same in every probe; probe.py's sampling and weights run unchanged
on it, with a stub in place of the solver. beta near ln(7.5) = 2.0 matches its growth."""
import argparse, math, os, random, re, shutil, sys, tempfile, threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe  # noqa: E402


def kids(nid, s):
    r = random.Random(nid)
    return [] if s <= 0 else [(f"{nid}.{i}", s - r.uniform(0.4, 1.6)) for i in range(r.randint(3, 12))]


def cpu_of(nid, s):
    return 10 + 50 * max(0.0, s) * random.Random(nid + "c").random()


def truth(nid, s):
    n, c = 1, cpu_of(nid, s)
    for k, ks in kids(nid, s):
        a, b = truth(k, ks)
        n, c = n + a, c + b
    return n, c


class Stub:
    """probe.Runner's interface, the solver replaced by the synthetic tree"""
    def __init__(self, k, beta, seed, run):
        self.a = argparse.Namespace(k=k, max_levels=40)
        self.cfg = {"k": k, "beta": beta, "seed": seed, "settings_file": os.path.abspath(__file__),
                    "commit": {"commit": "synthetic."}}
        self.run, self.map, self.stop = run, {}, threading.Event()

    def log(self, m):
        pass

    def task(self, wd, level, F=None, src=None, index=None):
        if F:
            nid, s = "r", float(F.split(":")[1])
        else:
            pnid, ps = self.map[os.path.dirname(src)]
            nid, s = kids(pnid, ps)[index]
        self.map[wd] = (nid, s)
        cpu = cpu_of(nid, s)
        tr = [{"cutoff": 100.0, "lb": 100 - s, "short": s, "iters": 300, "cut": 0, "cpu": cpu}]
        return {"children": [{"pred": random.Random(k).random()} for k, _ in kids(nid, s)], "trace": tr,
                "dnn_cpu": cpu, "task_cpu": cpu}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beta", type=float, default=1.8)
    ap.add_argument("--k", type=int, default=12)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--s0", type=float, default=2.6)
    a = ap.parse_args()
    tn, tc = truth("r", a.s0)
    est, cover, bounded = [], 0, 0
    for seed in range(a.reps):
        d = tempfile.mkdtemp()
        try:
            p = probe.Probe(Stub(a.k, a.beta, seed, d), "t", f"0:{a.s0}")
            p.go()
            line = open(os.path.join(p.dir, "summary.txt")).readline()
            c = float(re.search(r"DNN CPU ([\d,.]+) s", line).group(1).replace(",", ""))
            se = float(re.search(r"SE ([\d,.]+) s", line).group(1).replace(",", "")) if "SE" in line else 0.0
            est.append(c)
            cover += abs(c - tc) <= 1.96 * se
            bounded = len(p.nodes)
        finally:
            shutil.rmtree(d)
    m = sum(est) / len(est)
    sd = math.sqrt(sum((x - m) ** 2 for x in est) / (len(est) - 1))
    print(f"tree: {tn:,} nodes, {tc:,.0f} CPU s; {a.reps} probes, k {a.k}, beta {a.beta}, {bounded} nodes bounded each")
    print(f"mean {m:,.0f} = {m / tc:.3f} of the truth (SE of the mean {sd / math.sqrt(len(est)) / tc:.3f}); "
          f"one probe's sd {sd / tc:.2f} of the truth; +-1.96 SE held the truth {cover} of {a.reps} times")


if __name__ == "__main__":
    main()
