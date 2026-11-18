#!/usr/bin/env python3
"""A dives run end to end on the local backend (coord/dives.py, latest/estimate.C's survivor dives).

  python3 tests/dives_test.py [--instance nug17 --opt 1732] [--seconds 6]

1. nug17 (PARAM/default.par, the optimum as the incumbent, so the tree is fixed) solved serially:
   the true node count.
2. The same through the coordinator: 0.3-second units, 2 slots, for a few seconds; then paused,
   and its units left to return and merge.
3. A dives run from its queued units (3 seconds of survivor dives a unit), its coordinator and
   worker to the end, and `qapctl estimate`: the estimated tree must be within 15% of the count.
4. The paused run resumed to the end: its totals must equal the serial run's, so the dives run
   left its state files and parent states as they were.
5. With --flow: the same with nothing done by hand, `qapctl estimate-flow` (the solve for --seconds
   of its wall time, the pause, the dives run) and one `qapctl worker --follow`: both must exit 0
   with the solve paused, its dives run finished and the estimate within 15%.
A few minutes on an idle machine.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from coord.backend import open_backend  # noqa: E402


def qapctl(*args, **kw):
    return subprocess.run([sys.executable, "-m", "coord.qapctl"] + list(args), cwd=ROOT, check=True,
                          capture_output=True, text=True, **kw)


PROCS = []


def popen(*args, log):
    p = subprocess.Popen([sys.executable, "-m", "coord.qapctl"] + list(args), cwd=ROOT,
                         stdout=open(log, "w"), stderr=subprocess.STDOUT)
    PROCS.append(p)
    return p


def wait_for(cond, timeout, what):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if cond():
            return
        time.sleep(0.5)
    raise SystemExit(f"FAIL: timed out waiting for {what}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instance", default="nug17")
    ap.add_argument("--opt", type=float, default=1732)
    ap.add_argument("--seconds", type=float, default=6)
    ap.add_argument("--flow", action="store_true", help="also qapctl estimate-flow with a worker --follow")
    ap.add_argument("--dive-b", default="0.4", help="the dives' mode: survivor b, or kQ for Knuth's dives")
    ap.add_argument("--dive-sample", default="1", help="the share of the open nodes sampled")
    ap.add_argument("--dive-weights", default="unif", choices=["unif", "qpb"])
    a = ap.parse_args()
    dat = os.path.join(ROOT, "qaplib", f"{a.instance}.dat")
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    full = subprocess.run([os.path.join(ROOT, "latest", "bnbdrv"), "-p", "PARAM/default.par", "-u", str(a.opt), dat],
                          cwd=os.path.join(ROOT, "latest"), capture_output=True, text=True, check=True, env=env).stdout
    true_nodes = int(re.search(r"Nodes:\s*(\d+)", full).group(1))
    print(f"{a.instance}: {true_nodes:,} nodes serially")
    failed = 0
    with tempfile.TemporaryDirectory() as d:
      try:
          be = f"local:{d}/store"
          b = open_backend(be)
          run, dv = "dt", "dt-dives"
          qapctl("init", run, "--backend", be, "--instance", dat, "--params", "PARAM/default.par", "--ub", str(a.opt),
                 "--slots", "2", "--budget", "0.3", "--target-unit-cpu", "0.3", "--frontier-depth", "1",
                 "--visibility", "private", env=env)
          w = popen("worker", run, "--backend", be, "--slots", "2", "--workdir", f"{d}/w1", log=f"{d}/w1.log")
          c = popen("coordinator", run, "--backend", be, log=f"{d}/c1.log")
          time.sleep(a.seconds)
          qapctl("pause", run, "--backend", be)
          busy = lambda: {k: v for k, v in b.count_by_status(run).items() if k in ("leased", "done", "returned") and v}
          wait_for(lambda: not busy(), 120, "the paused run's units to return and merge")
          for p in (c, w):
              p.terminate()
              p.wait(timeout=60)
          counts = b.count_by_status(run)
          queued = counts.get("queued", 0)
          if run_status(b, run) != "paused" or not queued:
              print(f"FAIL: the solve finished or left nothing queued ({counts}); use fewer --seconds")
              return 1
          print(f"paused: {counts}")
          qapctl("dives", run, "--as", dv, "--backend", be, "--budget", "3", "--dive-b", a.dive_b,
               "--dive-sample", a.dive_sample, "--dive-weights", a.dive_weights, env=env)
          w = popen("worker", dv, "--backend", be, "--slots", "2", "--workdir", f"{d}/w2", "--exit-when-idle", "5",
                    log=f"{d}/w2.log")
          c = popen("coordinator", dv, "--backend", be, log=f"{d}/c2.log")
          wait_for(lambda: run_status(b, dv) == "finished", 600, "the dives run to finish")
          w.wait(timeout=120)
          c.wait(timeout=60)
          out = qapctl("estimate", dv, "--backend", be, "--out", f"{d}/est").stdout
          e = json.load(open(f"{d}/est/estimate.json"))
          ratio = e["total"]["nodes"] / true_nodes
          print(out.strip().splitlines()[2] if len(out.strip().splitlines()) > 2 else out)
          print(f"estimate {e['total']['nodes']:,.0f} nodes against {true_nodes:,}: {ratio:.3f} "
                f"({e['dives']} dives in {e['dived']} of {e['strata']} open nodes)")
          if not 0.85 <= ratio <= 1.15:
              print("FAIL: the estimate is more than 15% off")
              failed += 1
          if e["strata_with_survivors_not_dived"]:
              print(f"note: {e['strata_with_survivors_not_dived']} open nodes with survivors got no dive")
          # the solve, resumed, must finish with the serial totals: the dives touched nothing of it
          qapctl("resume", run, "--backend", be)
          w = popen("worker", run, "--backend", be, "--slots", "2", "--workdir", f"{d}/w3", "--exit-when-idle", "20",
                    log=f"{d}/w3.log")
          c = popen("coordinator", run, "--backend", be, log=f"{d}/c3.log")
          wait_for(lambda: run_status(b, run) == "finished", 600, "the resumed solve to finish")
          w.wait(timeout=120)
          c.wait(timeout=60)
          res = b.get_bytes(f"runs/{run}/results.txt").decode()
          got = int(re.search(r"Nodes:\s*(\d+)", res).group(1))
          print(f"resumed solve: {got:,} nodes")
          if got != true_nodes:
              print(f"FAIL: the resumed solve has {got} nodes, the serial run {true_nodes}")
              failed += 1
      finally:
        for p in PROCS:
            if p.poll() is None:
                p.kill()
                p.wait()
    if a.flow:
        failed += flow(a, dat, env, true_nodes)
    print("PASS" if not failed else f"FAIL: {failed} checks")
    return 1 if failed else 0


def flow(a, dat, env, true_nodes):
    """estimate-flow and worker --follow: the same estimate with nothing done by hand"""
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        try:
            be = f"local:{d}/store"
            b = open_backend(be)
            run = "ft"
            qapctl("init", run, "--backend", be, "--instance", dat, "--params", "PARAM/default.par", "--ub", str(a.opt),
                   "--slots", "2", "--budget", "0.3", "--target-unit-cpu", "0.3", "--frontier-depth", "1",
                   "--visibility", "private", env=env)
            w = popen("worker", run, "--backend", be, "--slots", "2", "--workdir", f"{d}/w", "--follow",
                      log=f"{d}/w.log")
            c = popen("estimate-flow", run, "--backend", be, "--solve-wall", str(a.seconds), "--dive-budget", "3",
                      "--dive-b", a.dive_b, "--dive-sample", a.dive_sample, "--dive-weights", a.dive_weights,
                      log=f"{d}/c.log")
            rc = c.wait(timeout=900)
            wrc = w.wait(timeout=300)
            st, sd = run_status(b, run), run_status(b, run + "-dives")
            print(f"flow: estimate-flow exit {rc}, worker exit {wrc}; {run} {st}, {run}-dives {sd}")
            if rc != 0 or wrc != 0 or st != "paused" or sd != "finished":
                print("FAIL: the flow did not end with the solve paused and its dives run finished")
                bad += 1
            else:
                e = json.loads(b.get_bytes(f"runs/{run}-dives/estimate.json"))
                ratio = e["total"]["nodes"] / true_nodes
                print(f"flow: estimate {e['total']['nodes']:,.0f} nodes against {true_nodes:,}: {ratio:.3f}")
                if not 0.85 <= ratio <= 1.15:
                    print("FAIL: the flow's estimate is more than 15% off")
                    bad += 1
        finally:
            for p in PROCS:
                if p.poll() is None:
                    p.kill()
                    p.wait()
    return bad


def run_status(b, run):
    r = b.get_run(run)
    return r.get("status") if r else None


if __name__ == "__main__":
    sys.exit(main())
