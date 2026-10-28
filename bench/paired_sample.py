#!/usr/bin/env python3
"""Compare params files on the same random sample of subtrees.

  python3 bench/paired_sample.py qaplib/nug25.dat --ub 3744 --depth 4 --samples 300 \\
      PARAM/default.par /path/to/variant.par [...]

A frontier pass with the first params file solves the tree down to depth D
(-d D -o F). The same random sample of its open nodes is then solved to the
end under every params file (-r, upper bound fixed). Because each variant
solves the same subtrees, most of the subtree-to-subtree variation cancels:
the ratio of summed CPU (variant / first) is far tighter than comparing two
independent sampling estimates, and is reported with a 95% bootstrap
interval.

This measures what a params change does below depth D. Changes that act at
depth <= D (the frontier itself) are not measured; the frontier is the
first file's. For nug24/25 the levels below depth 4 carry ~97% of the time.
Runs --jobs subtrees at a time (default 1); time ratios are only fair if
every variant sees the same load, so variants of one node run back to back.
"""

import argparse
import concurrent.futures as cf
import os
import random
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bench"))
from parallel import State  # noqa: E402

SOLVER = os.path.join(ROOT, "latest", "bnbdrv")


def totals(text):
    t = re.search(r"Total Time:\s*([\d.eE+-]+)", text)
    m = re.search(r"Nodes:\s*(\d+)", text)
    f = re.search(r"fw iterations:\s*(\d+)", text)
    if not (t and m):
        return None
    return int(m.group(1)), float(t.group(1)), int(f.group(1)) if f else 0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("instance")
    p.add_argument("params", nargs="+", help="params files (relative to latest/ or absolute); "
                   "the first is the baseline and builds the frontier")
    p.add_argument("--ub", required=True)
    p.add_argument("--depth", type=int, default=4)
    p.add_argument("--samples", type=int, default=300)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--solver", default=SOLVER)
    p.add_argument("--frontier", help="reuse this frontier state file")
    p.add_argument("--workdir")
    args = p.parse_args()

    solver = os.path.abspath(args.solver)
    sdir = os.path.dirname(solver)
    inst = os.path.abspath(args.instance)
    work = args.workdir or tempfile.mkdtemp(prefix="qap-paired-")
    os.makedirs(work, exist_ok=True)

    frontier = args.frontier or os.path.join(work, "frontier.state")
    if not args.frontier:
        r = subprocess.run([solver, "-p", args.params[0], "-u", str(args.ub),
                            "-d", str(args.depth), "-o", frontier, inst],
                           cwd=sdir, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"frontier pass failed:\n{r.stdout[-2000:]}{r.stderr}")
    st = State(frontier)
    n_open = len(st.nodes)
    print(f"frontier: depth {args.depth}, {st.counters[0]} nodes, {st.times[0]:.1f} s CPU, "
          f"{n_open} open nodes", flush=True)
    picks = random.Random(args.seed).sample(range(n_open), min(args.samples, n_open))
    template = State(frontier)
    template.zero_stats()

    def solve_node(k):
        """All params files on sampled node k, back to back."""
        i = picks[k]
        task = os.path.join(work, f"node{k}.state")
        template.write(task, [st.nodes[i]])
        res = []
        for par in args.params:
            out = subprocess.run([solver, "-p", par, "-u", str(args.ub), "-r", task, inst],
                                 cwd=sdir, capture_output=True, text=True).stdout
            t = totals(out)
            if t is None:
                raise RuntimeError(f"node {i} with {par} failed:\n{out[-1500:]}")
            res.append(t)
        os.remove(task)
        return k, res

    rows = [None] * len(picks)
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for done, fut in enumerate(cf.as_completed([ex.submit(solve_node, k)
                                                    for k in range(len(picks))]), 1):
            k, res = fut.result()
            rows[k] = res
            if done % 25 == 0 or done == len(picks):
                print(f"  {done}/{len(picks)} subtrees", flush=True)

    with open(os.path.join(work, "paired.csv"), "w") as f:
        f.write("node," + ",".join(f"nodes{j},cpu{j},fw{j}" for j in range(len(args.params))) + "\n")
        for k, res in enumerate(rows):
            f.write(f"{picks[k]}," + ",".join(f"{n},{c},{w}" for n, c, w in res) + "\n")

    rng = random.Random(args.seed + 1)
    base_n = sum(r[0][0] for r in rows)
    base_c = sum(r[0][1] for r in rows)
    print(f"\n{len(rows)} subtrees of {n_open}; baseline {args.params[0]}: "
          f"{base_n:,} nodes, {base_c:.1f} s")
    print(f"{'params':40s} {'nodes':>8s} {'CPU':>8s}   95% interval (CPU)")
    for j, par in enumerate(args.params):
        rn = sum(r[j][0] for r in rows) / base_n
        rc = sum(r[j][1] for r in rows) / base_c
        boots = []
        for _ in range(2000):
            idx = [rng.randrange(len(rows)) for _ in rows]
            b = sum(rows[i][0][1] for i in idx)
            boots.append(sum(rows[i][j][1] for i in idx) / b if b > 0 else 1.0)
        boots.sort()
        print(f"{os.path.basename(par):40s} {rn:8.3f} {rc:8.3f}   "
              f"[{boots[50]:.3f}, {boots[1949]:.3f}]")
    print(f"\nper-node results: {os.path.join(work, 'paired.csv')}")


if __name__ == "__main__":
    main()
