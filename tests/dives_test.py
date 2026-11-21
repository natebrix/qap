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
   of its wall time, the pause, the dives run; its decision point off) and one `qapctl worker
   --follow`: both must exit 0 with the solve paused, its dives run finished and the estimate within
   15%. Then (6) the decision point.
6. With --flow or --decide: estimate-flow with the decision point at 1% (`--decide 0.01`: an
   estimate of a tree this small always costs more), checked every second, 8 s of dives a unit:
   once the pilot's units (half of so few, coord/dives.py pilot_size) have merged the dives must stop and the solve go on to the end on the same worker, which comes
   back from the dives run: both exit 0, the solve finished with the serial run's nodes, the dives
   run stopped with the decision, and the readout saying so.
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
    ap.add_argument("--flow", action="store_true", help="also qapctl estimate-flow with a worker --follow (5, 6)")
    ap.add_argument("--decide", action="store_true", help="only the decision point's flow (6)")
    ap.add_argument("--dive-b", default="0.4", help="the dives' mode: survivor b, or kQ for Knuth's dives")
    ap.add_argument("--dive-sample", default="1", help="the share of the open nodes sampled")
    ap.add_argument("--dive-weights", default="unif", choices=["unif", "qpb"])
    ap.add_argument("--dive-env", action="append", default=[], help="QAP_EST_X=V for every dives unit (repeatable)")
    a = ap.parse_args()
    dat = os.path.join(ROOT, "qaplib", f"{a.instance}.dat")
    env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    full = subprocess.run([os.path.join(ROOT, "latest", "bnbdrv"), "-p", "PARAM/default.par", "-u", str(a.opt), dat],
                          cwd=os.path.join(ROOT, "latest"), capture_output=True, text=True, check=True, env=env).stdout
    true_nodes = int(re.search(r"Nodes:\s*(\d+)", full).group(1))
    print(f"{a.instance}: {true_nodes:,} nodes serially")
    failed = 0
    if a.decide:
        failed = decide_flow(a, dat, env, true_nodes)
        print("PASS" if not failed else f"FAIL: {failed} checks")
        return 1 if failed else 0
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
               "--dive-sample", a.dive_sample, "--dive-weights", a.dive_weights,
               *[x for kv in a.dive_env for x in ("--dive-env", kv)], env=env)
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
        failed += decide_flow(a, dat, env, true_nodes)
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
                      "--decide", "0", *[x for kv in a.dive_env for x in ("--dive-env", kv)], log=f"{d}/c.log")
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
                ro = (b.get_bytes(f"runs/{run}/readout.md") or b"").decode()
                print("flow: readout: " + (ro.splitlines()[2] if len(ro.splitlines()) > 2 else "missing"))
                if "**COMPLETE**: the estimate" not in ro or "## The estimate" not in ro:
                    print("FAIL: the source's readout does not give the estimate")
                    bad += 1
                if not 0.85 <= ratio <= 1.15:
                    print("FAIL: the flow's estimate is more than 15% off")
                    bad += 1
                # the paused solve's RESULTS block (coordinator.write_partial), its nodes the solve's own
                part = (b.get_bytes(f"runs/{run}/results_partial.txt") or b"").decode()
                st_nodes = json.loads(b.get_bytes(f"runs/{run}/status.json") or b"{}").get("nodes")
                m = re.search(r"Nodes: (\d+) /", part)
                print("flow: results_partial.txt: " + (part.splitlines()[0][:100] if part else "missing")
                      + (f"; Nodes {m.group(1)}, status.json {st_nodes}" if m else ""))
                if not part.startswith("PARTIAL: the solve paused") or "---- RESULTS ----" not in part or \
                        "\nlevel " not in part or not m or int(m.group(1)) != st_nodes:
                    print("FAIL: the paused solve's results_partial.txt is missing, malformed or not its totals")
                    bad += 1
        finally:
            for p in PROCS:
                if p.poll() is None:
                    p.kill()
                    p.wait()
    return bad


def decide_flow(a, dat, env, true_nodes):
    """the decision point (coord/dives.py): the dives stop for the solve, which finishes on the same worker"""
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        try:
            be = f"local:{d}/store"
            b = open_backend(be)
            run = "dc"
            qapctl("init", run, "--backend", be, "--instance", dat, "--params", "PARAM/default.par", "--ub", str(a.opt),
                   "--slots", "2", "--budget", "0.3", "--target-unit-cpu", "0.3", "--frontier-depth", "1",
                   "--visibility", "private", env=env)
            w = popen("worker", run, "--backend", be, "--slots", "2", "--workdir", f"{d}/w", "--follow",
                      log=f"{d}/w.log")
            # dives of 8 s a unit, so that the units merge a few at a time and a check falls between
            # (with 3 s all of them can merge in one tick, and the run finishes before any check)
            c = popen("estimate-flow", run, "--backend", be, "--solve-wall", str(a.seconds), "--dive-budget", "8",
                      "--dive-b", a.dive_b, "--dive-sample", a.dive_sample, "--dive-weights", a.dive_weights,
                      "--decide", "0.01", "--decide-every", "1",
                      *[x for kv in a.dive_env for x in ("--dive-env", kv)], log=f"{d}/c.log")
            rc = c.wait(timeout=900)
            wrc = w.wait(timeout=300)
            st, sd = run_status(b, run), run_status(b, run + "-dives")
            dv = b.get_run(run + "-dives") or {}
            dec = dv.get("decision") or {}
            print(f"decide: estimate-flow exit {rc}, worker exit {wrc}; {run} {st}, {run}-dives {sd}, "
                  f"decision {dec.get('verdict')} after {dec.get('merged')} of {dec.get('units')} units, "
                  f"pilot {dec.get('pilot_merged')} of {dec.get('pilot')} (share {dec.get('share')})")
            if rc != 0 or wrc != 0 or st != "finished" or sd != "stopped" or dec.get("verdict") != "solve":
                print("FAIL: the flow did not stop the dives for the solve and finish it")
                bad += 1
            else:
                res = b.get_bytes(f"runs/{run}/results.txt").decode()
                got = int(re.search(r"Nodes:\s*(\d+)", res).group(1))
                print(f"decide: the solve: {got:,} nodes")
                if got != true_nodes:
                    print(f"FAIL: the solve after the decision has {got} nodes, the serial run {true_nodes}")
                    bad += 1
                if not dec.get("pilot") or dec.get("pilot_merged") != dec.get("pilot"):
                    print(f"FAIL: it decided before its pilot had merged ({dec.get('pilot_merged')} of {dec.get('pilot')})")
                    bad += 1
                ro = (b.get_bytes(f"runs/{run}/readout.md") or b"").decode()
                print("decide: readout: " + (ro.splitlines()[2] if len(ro.splitlines()) > 2 else "missing"))
                if "**COMPLETE**: solved" not in ro or "DECIDED TO SOLVE" not in ro:
                    print("FAIL: the readout does not say the estimate decided to solve")
                    bad += 1
                if "back to" not in open(f"{d}/w.log").read():
                    print("FAIL: the worker did not come back from the dives run")
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
