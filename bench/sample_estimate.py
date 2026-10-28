#!/usr/bin/env python3
"""Estimate the size and CPU time of a full solve by sampling subtrees.

  python3 bench/sample_estimate.py qaplib/nug25.dat --ub 3744 --depth 3 --samples 40

1. A frontier pass solves the tree down to depth D (-d D -o F), like
   bench/parallel.py. Its nodes and CPU time are counted exactly.
2. A random sample of the open nodes at depth D is solved to the end, each
   on its own (-r), with the upper bound fixed (use the optimum, so the
   work below a node does not depend on the order).
3. Total = frontier + (open nodes) x (sample mean), for nodes and CPU
   seconds, with a 95% interval from the sample's standard error. Subtree
   sizes are heavy-tailed, so the interval is rough; more samples or a
   deeper frontier narrow it.

Validate on an instance you can also solve fully before trusting it on one
you can't.
"""

import argparse
import math
import os
import random
import re
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
from parallel import State  # noqa: E402

SOLVER = os.path.join(ROOT, "latest", "bnbdrv")


def totals(text):
    t = re.search(r"Total Time:\s*([\d.eE+-]+)", text)
    m = re.search(r"Nodes:\s*(\d+)", text)
    return (int(m.group(1)) if m else None, float(t.group(1)) if t else None)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("instance")
    p.add_argument("--ub", required=True, help="upper bound (use the optimum)")
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--samples", type=int, default=40)
    p.add_argument("--params", default="PARAM/default.par", help="relative to latest/")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--solver", default=SOLVER)
    p.add_argument("--workdir", help="keep files here (default: a temp dir)")
    args = p.parse_args()

    solver = os.path.abspath(args.solver)
    sdir = os.path.dirname(solver)
    inst = os.path.abspath(args.instance)
    work = args.workdir or tempfile.mkdtemp(prefix="qap-est-")
    os.makedirs(work, exist_ok=True)
    base = [solver, "-p", args.params, "-u", str(args.ub)]

    frontier = os.path.join(work, "frontier.state")
    t0 = time.time()
    r = subprocess.run(base + ["-d", str(args.depth), "-o", frontier, inst],
                       cwd=sdir, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"frontier pass failed (rc={r.returncode}):\n{r.stdout[-2000:]}{r.stderr}")
    st = State(frontier)
    f_nodes, f_cpu = st.counters[0], st.times[0]
    n_open = len(st.nodes)
    print(f"frontier: depth {args.depth}, {f_nodes} nodes, {f_cpu:.1f} s CPU, "
          f"{n_open} open nodes ({time.time() - t0:.0f} s wall)", flush=True)
    if n_open == 0:
        print(f"estimate: {f_nodes} nodes, {f_cpu:.1f} s (solved in the frontier pass)")
        return

    rng = random.Random(args.seed)
    picks = rng.sample(range(n_open), min(args.samples, n_open))
    template = State(frontier)
    template.zero_stats()
    nodes, cpus = [], []
    for k, i in enumerate(picks):
        task = os.path.join(work, f"sample{k}.state")
        template.write(task, [st.nodes[i]])
        out = subprocess.run(base + ["-r", task, inst], cwd=sdir,
                             capture_output=True, text=True).stdout
        n, c = totals(out)
        if n is None:
            sys.exit(f"sample {k} failed:\n{out[-2000:]}")
        nodes.append(n)
        cpus.append(c)
        os.remove(task)
        print(f"  sample {k + 1}/{len(picks)}: node {i}: {n} nodes, {c:.2f} s", flush=True)

    def est(xs, base_value):
        k = len(xs)
        mean = sum(xs) / k
        sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (k - 1)) if k > 1 else 0.0
        fpc = math.sqrt((n_open - k) / (n_open - 1)) if n_open > 1 else 0.0
        half = 1.96 * sd / math.sqrt(k) * fpc * n_open
        return base_value + n_open * mean, half

    en, hn = est(nodes, f_nodes)
    ec, hc = est(cpus, f_cpu)
    print(f"\nestimate from {len(picks)} of {n_open} subtrees at depth {args.depth}:")
    print(f"  nodes {en:,.0f}  (95% interval +-{hn:,.0f})")
    print(f"  CPU   {ec:,.0f} s = {ec / 3600:.2f} h  (95% interval +-{hc:,.0f} s)")
    print(f"  largest sampled subtree: {max(nodes):,} nodes, {max(cpus):.1f} s")


if __name__ == "__main__":
    main()
