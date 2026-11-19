#!/usr/bin/env python3
"""Validity runs of dom.pair through the solver (README, "Validity").

For each instance: the reference optimum (QAPSOLN, else brute force for n <= 11,
else the key-off run), then bnbdrv with dom.pair = 0 and = 1 from the same
incumbent, each to the end (-o: the final state must have 0 open nodes).
Incumbent modes:
  opt1   the optimum + 1, rounded up by the solver to the objective grid, so
         the search must find the optimum itself
  poor   2 x optimum + 100 (far above: a poor starting solution)
  none   no incumbent at all
One JSON line per run to the output file. One solver process at a time.

  python3 run_valid.py OUT.jsonl MODE TIMEOUT inst [inst ...]
"""
import json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
WT = "/home/user/wt-dombuild"
SOLVER = os.environ.get("DOM_SOLVER", f"{WT}/latest/bnbdrv")
QAPLIB = f"{WT}/qaplib"
SCR = os.environ.get("DOM_SCRATCH", "/tmp/claude-0/-home-user/92b66c47-7a9f-5ca5-8a65-57e14ac5b075/scratchpad")
PAR = {0: f"{HERE}/par/default.par", 1: f"{HERE}/par/default_dom.par"}
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")


def ref_opt(name, n):
    s = f"{WT}/QAPSOLN/{name}.sln"
    if os.path.exists(s):
        return int(open(s).read().split()[1]), "QAPSOLN", None
    t = open(f"{QAPLIB}/{name}.dat").read().split()
    if len(t) - 1 - 2 * n * n not in (0, n * n):
        return None, None, None   # trailing tokens the solver reads its own way (lipa10*): the key-off run decides
    if n <= 11:
        out = subprocess.run([f"{HERE}/bin/bf", f"{QAPLIB}/{name}.dat"], capture_output=True, text=True).stdout.split()
        return int(out[0]), "bruteforce", int(out[1])
    return None, None, None


def run(name, key, ub, timeout, trace=None):
    st = f"{SCR}/{name}.{key}.state"
    cmd = [SOLVER, "-p", PAR[key], "-o", st]
    if ub is not None:
        cmd += ["-u", str(ub)]
    cmd.append(f"{QAPLIB}/{name}.dat")
    env = dict(ENV)
    if trace:
        if os.path.exists(trace): os.remove(trace)
        env["QAP_NODE_TRACE"] = trace
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=f"{WT}/latest", capture_output=True, text=True, timeout=timeout, env=env)
        text, rc, to = r.stdout + r.stderr, r.returncode, False
    except subprocess.TimeoutExpired as e:
        text = (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        rc, to = None, True
    res = {"instance": name, "key": key, "ub": ub, "wall": round(time.time() - t0, 2)}
    m = re.search(r"Nodes:\s*(\d+)\s*/\s*(\d+)\s*/\s*(\d+)", text)
    if m: res["nodes"], res["fathom"], res["enumerate"] = map(int, m.groups())
    found = re.findall(r"New objective value:\s*([0-9.eE+-]+)", text)
    res["found"] = float(found[-1]) if found else None
    m = re.search(r"Pair exchange dominance: (\d+) nodes dropped", text)
    res["dropped"] = int(m.group(1)) if m else (0 if key == 0 else None)
    m = re.search(r"dom.pair: the test is off for this instance: (.*)", text)
    res["dom_unavailable"] = m.group(1) if m else None
    res["asym_warning"] = "neither A nor B is symmetric" in text
    res["sym_file"] = "read symmetry info" in text
    m = re.search(r"incumbent \S+ is not an attainable objective value; using (\S+)", text)
    res["ub_used"] = float(m.group(1)) if m else ub
    res["finished"] = "---- RESULTS ----" in text
    res["status"] = "timeout" if to else ("ok" if rc == 0 and res["finished"] else f"error(rc={rc})")
    if os.path.exists(st) and not to:
        hdr = open(st).read(4000).split("\n")
        res["state_version"] = hdr[0]
        inc = [l for l in hdr if l.startswith("incumbent ")]
        res["final_incumbent"] = float(inc[0].split()[1]) if inc else None
        with open(st) as f:
            for l in f:
                if l.startswith("nodes "):
                    res["open_nodes"] = int(l.split()[1])
        os.remove(st)
    return res


def main():
    out, mode, timeout = sys.argv[1], sys.argv[2], float(sys.argv[3])
    for name in sys.argv[4:]:
        n = int(open(f"{QAPLIB}/{name}.dat").read().split()[0])
        opt, src, nopt = ref_opt(name, n)
        if opt is None:   # the key-off run without an incumbent is the reference
            r0 = run(name, 0, None, timeout)
            opt, src = (int(r0["final_incumbent"]), "solver(off,no ub)") if r0["status"] == "ok" else (None, None)
        if opt is None:
            print(name, "no reference optimum", flush=True); continue
        ub = {"opt1": opt + 1, "poor": 2 * abs(opt) + 100, "none": None}[mode]
        for key in (0, 1):
            r = run(name, key, ub, timeout)
            r.update(n=n, opt=opt, opt_src=src, n_optima=nopt, mode=mode)
            ok = (r["status"] == "ok" and r.get("open_nodes") == 0 and r.get("final_incumbent") == opt)
            r["valid"] = ok
            with open(out, "a") as f: f.write(json.dumps(r) + "\n")
            print(f"{name:10s} key={key} {r['status']:8s} nodes={r.get('nodes')} dropped={r.get('dropped')} "
                  f"final={r.get('final_incumbent')} opt={opt} open={r.get('open_nodes')} "
                  f"{'OK' if ok else 'CHECK'} {r['wall']}s {r.get('dom_unavailable') or ''}", flush=True)


if __name__ == "__main__":
    main()
