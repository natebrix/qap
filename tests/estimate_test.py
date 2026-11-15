#!/usr/bin/env python3
"""Check the tree-size estimator (bnbdrv -e, latest/estimate.C; thesis chapter 4).

  python3 tests/estimate_test.py [--solver PATH]

nug15 is solved once to count its tree. Then the estimator runs 3 x 400 dives
each way a dive can weight children: by their predicted bounds (the default)
and by their own QPB bounds (QAP_EST_CHILD_QPB=1). Knuth's estimate is unbiased
whatever the weights, so both pooled means must come within 15% of the count,
and so must their CPU estimates of the run's CPU (a weighting that changed the
search itself, or counted its own bounds as the tree's time, fails here: the
first try at QAP_EST_CHILD_QPB read 0.16 of nug20's nodes). With the setting the
first draw must say it scored its nodes. About 20 seconds.
"""

import argparse
import os
import re
import statistics
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCE, UB = "nug15", 1150
SEEDS, DIVES, TOLERANCE = (1, 2, 3), 400, 0.15


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", default=os.path.join(ROOT, "latest", "bnbdrv"))
    a = ap.parse_args()
    lat = os.path.join(ROOT, "latest")
    common = ["-p", os.path.join(lat, "PARAM", f"{INSTANCE}.par"), "-s", os.path.join(lat, "SYM", f"{INSTANCE}.sym"),
              "-u", str(UB), os.path.join(ROOT, "qaplib", f"{INSTANCE}.dat")]
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")
    full = subprocess.run([a.solver] + common, cwd=lat, env=env, capture_output=True, text=True, check=True).stdout
    nodes = int(re.search(r"Nodes:\s*(\d+)", full).group(1))
    cpu = float(re.search(r"Total Time:\s*(\S+)", full).group(1))
    print(f"{INSTANCE}: {nodes:,} nodes, {cpu:.2f} CPU s")
    failed = 0
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "estimate.par"), "w") as f:   # start level, dives, print every, print all every, q, ignore depth
            f.write(f"2 {DIVES} 1 100000 2 1000\n")
        for arm in ("0", "1"):
            n_est, t_est = [], []
            for seed in SEEDS:
                hist = os.path.join(d, f"dives_{arm}_{seed}.txt")
                out = subprocess.run([a.solver, "-e"] + common, cwd=d, capture_output=True, text=True, check=True,
                                     env=dict(env, QAP_EST_CHILD_QPB=arm, QAP_EST_SEED=str(seed), QAP_EST_HIST=hist)).stdout
                if arm == "1" and "QAP_EST_CHILD_QPB: the first draw's" not in out:
                    print(f"FAIL: QAP_EST_CHILD_QPB=1, seed {seed}: the first draw did not score its nodes")
                    failed += 1
                rows = [l.split() for l in open(hist) if l.strip()]
                if len(rows) != DIVES:
                    print(f"FAIL: QAP_EST_CHILD_QPB={arm}, seed {seed}: {len(rows)} dives, not {DIVES}")
                    failed += 1
                n_est += [float(r[1]) for r in rows]
                t_est += [float(r[2]) for r in rows]
            rn, rt = statistics.fmean(n_est) / nodes, statistics.fmean(t_est) / cpu
            ok = abs(rn - 1) <= TOLERANCE and abs(rt - 1) <= TOLERANCE
            failed += not ok
            print(f"{'ok  ' if ok else 'FAIL'} QAP_EST_CHILD_QPB={arm}: {len(n_est)} dives, nodes {rn:.3f}, CPU {rt:.3f} of the run's")
    print("estimator: all checks passed" if not failed else f"estimator: {failed} checks failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
