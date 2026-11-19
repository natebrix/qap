#!/usr/bin/env python3
"""dom.pair = 1 through the solver's machinery (README, "Repo checks"):
  resume    straight run vs -t slices / -r resumes (state version 10 round trip):
            the same Nodes line, objective, per-level counts and dropped nodes
  parallel  bench/parallel.py --jobs 2 vs the serial run: the same Nodes line and
            objective, and the dropped count summed from the workers' states
  trace     QAP_NODE_TRACE: lines = node counter + popped + dropped, the X lines
            per depth = dom_f_level, each X node's fixes once
  python3 keyon_mech.py [resume] [parallel] [trace]
"""
import collections, os, re, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
WT = "/home/user/wt-dombuild"
SOLVER = f"{WT}/latest/bnbdrv"
PAR = f"{HERE}/par/default_dom.par"
CASES = [("nug15", 1150), ("had14", 2724), ("tai12b", 39464925), ("chr12c", 11156), ("scr15", 51140), ("nug14", None)]
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")


def run(args, env=ENV):
    r = subprocess.run([SOLVER] + args, cwd=f"{WT}/latest", capture_output=True, text=True, timeout=3600, env=env)
    return r.returncode, r.stdout + r.stderr


def summ(text):
    g = lambda p: (re.search(p, text).group(1) if re.search(p, text) else None)
    m = re.search(r"\nlevel .*?\n(.*?)\n\n", text, re.S)
    return {"nodes": g(r"Nodes:\s*(\d+\s*/\s*\d+\s*/\s*\d+)"), "dropped": g(r"Pair exchange dominance: (\d+)"),
            "per-level": " ".join(l.split()[1] for l in m.group(1).splitlines()) if m else None}


def state_arrays(path):
    out = {}
    for l in open(path):
        t = l.split()
        if t and t[0].endswith("_level") or (t and t[0] in ("counters",)):
            out[t[0]] = [float(x) for x in t[1:]]
        if t and t[0] in ("warm_starts", "nodes"): break
    return out


def resume():
    tmp = tempfile.mkdtemp(dir=HERE + "/logs")
    fails = 0
    for name, ub in CASES:
        dat = f"{WT}/qaplib/{name}.dat"; base = ["-p", PAR] + (["-u", str(ub)] if ub else [])
        rc, full = run(base + [dat]); a = summ(full)
        st = f"{tmp}/{name}.state"
        rc, out = run(base + ["-o", st, "-t", "0.2", dat]); n = 1; objs = re.findall(r"New objective value:\s*(\S+)", out)
        ver = set([open(st).readline().strip()])
        while rc == 3:
            ver.add(open(st).readline().strip())
            rc, out = run(["-p", PAR, "-o", st, "-t", "0.2", "-r", st, dat]); n += 1
            objs += re.findall(r"New objective value:\s*(\S+)", out)
        b = summ(out)
        fo = re.findall(r"New objective value:\s*(\S+)", full)
        ok = rc == 0 and a == b and (objs[-1:] == fo[-1:])
        fails += not ok
        print(f"resume {'ok  ' if ok else 'FAIL'} {name}: {n} slices, {a} vs {b}, objective {fo[-1:]} vs {objs[-1:]}, states {sorted(ver)}", flush=True)
        os.remove(st)
    os.rmdir(tmp)
    return fails


def parallel():
    fails = 0
    for name, ub in CASES[:5]:
        dat = f"{WT}/qaplib/{name}.dat"
        rc, full = run(["-p", PAR, "-u", str(ub), dat]); a = summ(full)
        r = subprocess.run([sys.executable, f"{WT}/bench/parallel.py", dat, "--jobs", "2", "--params", PAR,
                            "--solver", SOLVER, "--ub", str(ub)], capture_output=True, text=True, env=ENV, timeout=3600)
        out = r.stdout + r.stderr; b = summ(out)
        ok = r.returncode == 0 and a["nodes"] == b["nodes"]
        fails += not ok
        print(f"parallel {'ok  ' if ok else 'FAIL'} {name}: serial {a['nodes']} dropped {a['dropped']}; parallel rc={r.returncode} {b['nodes']} dropped {b['dropped']}"
              + ("" if ok else "\n" + out[-1500:]), flush=True)
    return fails


def trace():
    fails = 0
    tmp = tempfile.mkdtemp(dir=HERE + "/logs")
    for name, ub in CASES[:5]:
        dat = f"{WT}/qaplib/{name}.dat"; tr = f"{tmp}/{name}.trace"; st = f"{tmp}/{name}.state"
        rc, out = run(["-p", PAR, "-u", str(ub), "-o", st, dat], env=dict(ENV, QAP_NODE_TRACE=tr))
        lines = [l.split() for l in open(tr)]
        arr = state_arrays(st)
        node = int(arr["counters"][0]); pop = int(sum(arr["pop_f_level"])); dom = int(sum(arr["dom_f_level"]))
        xs = collections.Counter(int(t[1]) for t in lines if t[7] == "X")
        keys = collections.Counter(t[2] for t in lines)
        ok = (len(lines) == node + pop + dom and all(xs[d] == int(arr["dom_f_level"][d]) for d in range(len(arr["dom_f_level"])))
              and max(keys.values()) == 1 and sum(xs.values()) == dom)
        fails += not ok
        print(f"trace {'ok  ' if ok else 'FAIL'} {name}: {len(lines)} lines = node {node} + popped {pop} + dropped {dom}; X lines {sum(xs.values())}; "
              f"outcomes {dict(collections.Counter(t[7] for t in lines))}", flush=True)
        os.remove(tr); os.remove(st)
    os.rmdir(tmp)
    return fails


if __name__ == "__main__":
    what = sys.argv[1:] or ["resume", "parallel", "trace"]
    f = sum({"resume": resume, "parallel": parallel, "trace": trace}[w]() for w in what)
    print("FAILURES", f)
    sys.exit(1 if f else 0)
