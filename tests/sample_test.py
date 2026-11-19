#!/usr/bin/env python3
"""A sampled A/B end to end on the local backend (coord/sample.py, experiments/estimates_nov19/paired.py).

  python3 tests/sample_test.py [--seconds 12]

1. nug24 under the nug30 recipe (the DNN at every depth, its states shipped: the sampled nodes'
   parents' states must come along) through the coordinator, 2-second units, traced, for --seconds;
   then paused, its units back.
2. A dives run of it (Knuth's dives, 2 seconds a unit) to the end: the size measure.
3. Two sampled runs, arms a and b with the same params, each of 4 open nodes at depth 2 and 4 at depth 3
   drawn by the dives run's estimates, each solved to the end. Checks: the two keep the same sample;
   every sampled node's subtree has the same nodes in both (the same params, the same parents' states:
   the same trees); the arms' traced nodes are all under the sampled nodes; the source is still paused
   with all its units queued (a sample reads its states, never changes them); paired.py reads both.
About two minutes on an idle machine.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from coord.backend import open_backend  # noqa: E402

PROCS = []
PAR = "../infra/phase1/recipes/nug30.par"


def qapctl(*args, env=None):
    r = subprocess.run([sys.executable, "-m", "coord.qapctl"] + list(args), cwd=ROOT, capture_output=True, text=True,
                       env=env)
    if r.returncode != 0:
        raise SystemExit(f"FAIL: qapctl {' '.join(args)}:\n{r.stdout[-3000:]}{r.stderr[-3000:]}")
    return r.stdout


def popen(*args, log, env=None):
    p = subprocess.Popen([sys.executable, "-m", "coord.qapctl"] + list(args), cwd=ROOT, stdout=open(log, "w"),
                         stderr=subprocess.STDOUT, env=env)
    PROCS.append(p)
    return p


def wait_for(cond, timeout, what):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if cond():
            return
        time.sleep(0.5)
    raise SystemExit(f"FAIL: timed out waiting for {what}")


def status(b, run):
    r = b.get_run(run)
    return r.get("status") if r else None


def run_to_end(b, be, run, d, env):
    w = popen("worker", run, "--backend", be, "--slots", "2", "--workdir", f"{d}/w-{run}", "--exit-when-idle", "10",
              log=f"{d}/w-{run}.log", env=env)
    c = popen("coordinator", run, "--backend", be, log=f"{d}/c-{run}.log", env=env)
    wait_for(lambda: status(b, run) == "finished", 900, f"{run} to finish")
    c.wait(timeout=60)
    w.wait(timeout=120)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seconds", type=float, default=12)
    a = ap.parse_args()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OMP_NUM_THREADS="1")
    dat = os.path.join(ROOT, "qaplib", "nug24.dat")
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        try:
            be = f"local:{d}/store"
            b = open_backend(be)
            unit = ["--slots", "2", "--budget", "2", "--target-unit-cpu", "2", "--trace", "--visibility", "private"]
            qapctl("init", "src", "--backend", be, "--instance", dat, "--params", PAR, "--ub", "3488",
                   "--frontier-depth", "1", *unit, env=env)
            w = popen("worker", "src", "--backend", be, "--slots", "2", "--workdir", f"{d}/w-src", log=f"{d}/w-src.log", env=env)
            c = popen("coordinator", "src", "--backend", be, log=f"{d}/c-src.log", env=env)
            time.sleep(a.seconds)
            qapctl("pause", "src", "--backend", be)
            busy = lambda: {k: v for k, v in b.count_by_status("src").items() if k in ("leased", "done", "returned") and v}
            wait_for(lambda: not busy(), 300, "the paused source's units to return")
            for p in (c, w):
                p.terminate()
                p.wait(timeout=60)
            queued = len(b.units_by_status("src", "queued"))
            if status(b, "src") != "paused" or not queued:
                print(f"FAIL: the source finished or left nothing queued ({b.count_by_status('src')}): fewer --seconds")
                return 1
            print(f"source paused: {b.count_by_status('src')}")
            qapctl("dives", "src", "--as", "src-dives", "--backend", be, "--budget", "2", "--dive-b", "k1.5", env=env)
            run_to_end(b, be, "src-dives", d, env)
            print("dives run finished")
            for arm in ("a", "b"):
                out = qapctl("init", arm, "--backend", be, "--instance", dat, "--params", PAR, "--ub", "3488",
                             "--sample-from", "src", "--sample", "2:4,3:4", *unit, env=env)
                print(f"{arm}: " + next((l for l in out.splitlines() if "sample:" in l), "").split("sample: ", 1)[-1])
                run_to_end(b, be, arm, d, env)
            sa, sb = (json.loads(b.get_bytes(f"runs/{r}/sample.json")) for r in ("a", "b"))
            if sa != sb:
                print("FAIL: the two arms keep different samples")
                bad += 1
            n = len(sa["nodes"])
            print(f"sample: {n} nodes, weights {', '.join(f'{x['w']:.1f}' for x in sa['nodes'])}")
            if n == 0:
                print("FAIL: an empty sample")
                bad += 1
            sys.path.insert(0, os.path.join(ROOT, "experiments", "estimates_nov19"))
            import paired
            roots = {x["key"]: x["depth"] for x in sa["nodes"]}
            acc = {}
            for r in ("a", "b"):
                _, _, lines, _ = paired.read_store(b, r)
                acc[r] = paired.subtrees(lines, roots)
                traced = sum(1 for l in lines if l.startswith("node "))
                under = sum(v[1] for v in acc[r].values())
                if traced != under:
                    print(f"FAIL: {r}: {traced} traced nodes, {under} under the sampled ones")
                    bad += 1
            diff = [k for k in roots if acc["a"][k][1] != acc["b"][k][1]]
            print(f"subtrees: " + ", ".join(f"{acc['a'][k][1]}" for k in roots) + f" nodes in a; {len(diff)} differ in b")
            if diff:
                print(f"FAIL: {len(diff)} sampled nodes' subtrees differ between two arms with the same params")
                bad += 1
            if status(b, "src") != "paused" or len(b.units_by_status("src", "queued")) != queued:
                print("FAIL: the source changed")
                bad += 1
            out = subprocess.run([sys.executable, os.path.join(ROOT, "experiments", "estimates_nov19", "paired.py"),
                                  "--backend", be, "a", "b", "--reps", "200"], cwd=ROOT, capture_output=True, text=True)
            print(out.stdout.strip())
            if out.returncode != 0 or "1.00x" not in out.stdout:
                print(f"FAIL: paired.py\n{out.stderr[-2000:]}")
                bad += 1
        finally:
            for p in PROCS:
                if p.poll() is None:
                    p.kill()
                    p.wait()
    print("PASS" if not bad else f"FAIL: {bad} checks")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
