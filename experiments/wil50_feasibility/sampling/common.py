"""Shared pieces of the sampling launchers (README.md): the full run's
configuration from fullrun/settings.sh and fullrun.par, one-node state files,
solver tasks pinned to one CPU, the DNN trace format and the seeded hash that
draws the stage-B sample."""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

D = os.path.dirname(os.path.abspath(__file__))
R = os.path.abspath(os.path.join(D, "../../.."))
sys.path.insert(0, os.path.join(R, "bench"))
from parallel import State, Warm  # noqa: E402,F401  (the solver's state-file format, latest/state.h)
import envpar  # noqa: E402  (settings overrides into a params file's [settings] section)

# SAMPLING_SETTINGS / SAMPLING_PARAMS: another instance's settings file and params file
# (experiments/estimates/INSTANCE/settings.sh, which defines the same FULLRUN_* names and
# fullrun_check); unset, the full wil50 run's
SETTINGS = os.environ.get("SAMPLING_SETTINGS") or os.path.join(R, "experiments", "wil50_feasibility", "fullrun", "settings.sh")
# the full run's params file (fullrun/README.md): nug_deep_xbar9_all.par (rule 9, lookahead
# between the two DNN-marginal lines, through depth 6, rule 2 below) with every solver
# setting in its [settings] section, the depth spec included (dnn.depths = 1:4, the DNN
# bound from depth 1, roadmap/tree_size.md section 4, through depth 4, then dnn.ext = 20
# levels below a miss: from depth 3 down exactly the absolute depths of the measured
# node-10 arms; dnn.giveup = 300; dnn.warm = 1). Until Nov 9 the settings were the
# environment of settings.sh and the depth spec the flag -x 1:4:1:300:1.
PARAMS = os.environ.get("SAMPLING_PARAMS") or os.path.join(R, "experiments", "wil50_feasibility", "fullrun", "fullrun.par")
XSPEC = ""   # no -x: the params file carries the depth spec (a value here overrides it)
SOLVER_PREFIXES = ("QAP_", "DNN_")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def log_line(path, msg):
    line = f"{time.strftime('%F %T')} {msg}"
    print(line, flush=True)
    if path:
        with open(path, "a") as f:
            f.write(line + "\n")


def fullrun_env(binary):
    """(env, info): exactly the solver variables fullrun/settings.sh exports (the
    thread variables and the diagnostics; the caller's own QAP_*/DNN_* variables
    are dropped), after its fullrun_check on `binary`; info holds
    FULLRUN_INSTANCE, FULLRUN_SYM, FULLRUN_UB and FULLRUN_PAR"""
    base = {k: v for k, v in os.environ.items() if not k.startswith(SOLVER_PREFIXES)}
    script = ('. "$1" || exit 1; fullrun_check "$2" >&2 || exit 1; '
              'printf "FULLRUN_INSTANCE=%s\\0FULLRUN_SYM=%s\\0FULLRUN_UB=%s\\0FULLRUN_PAR=%s\\0--\\0" '
              '"$FULLRUN_INSTANCE" "$FULLRUN_SYM" "$FULLRUN_UB" "$FULLRUN_PAR"; env -0')
    r = subprocess.run(["bash", "-c", script, "bash", SETTINGS, binary], env=base, capture_output=True)
    if r.returncode != 0:
        raise SystemExit("fullrun_check failed:\n" + r.stderr.decode())
    head, _, rest = r.stdout.decode().partition("--\0")
    info = dict(x.split("=", 1) for x in head.split("\0") if "=" in x)
    env = {}
    for x in rest.split("\0"):
        if "=" not in x:
            continue
        k, v = x.split("=", 1)
        if k.startswith(SOLVER_PREFIXES) or k in THREAD_VARS:
            env[k] = v
    return env, info


def params_with(base, overrides, out, binary=None):
    """`out`: the params file `base` with the settings among `overrides` (keys, or
    the pre-Phase-2 environment names) in its [settings] section; returns the
    overrides that are not settings (diagnostics, thread variables), for the
    environment. An existing `out` is kept (a resumed run keeps its file)."""
    if os.path.exists(out):
        settings, residual = envpar.split_env(overrides or {}, envpar.settings_map(binary or envpar.SOLVER))
        return residual
    return envpar.write_par(base, overrides or {}, out, solver=binary or envpar.SOLVER)


def dnn_flags(xspec):
    """the -x flag when a depth spec is given (the params file has it otherwise)"""
    return ["-x", xspec] if xspec else []


def solver_env(settings_env, overrides=None):
    """the process environment for a solver: the caller's minus any QAP_*/DNN_*,
    plus the settings, plus overrides (a value of None removes the variable)"""
    e = {k: v for k, v in os.environ.items() if not k.startswith(SOLVER_PREFIXES)}
    e.update(settings_env)
    for k, v in (overrides or {}).items():
        if v is None:
            e.pop(k, None)
        else:
            e[k] = v
    return e


def private_binary(src, bindir):
    """a private copy of the solver (a rebuild during the run must not change it)
    and the commit it was built from"""
    os.makedirs(bindir, exist_ok=True)
    dst = os.path.join(bindir, "bnbdrv")
    info = os.path.join(bindir, "commit.txt")
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=R, capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "latest"], cwd=R).returncode != 0
        newer = [f for f in ("dnn.C", "dnncut.C", "dnn.h", "bnb.C", "branch.C", "state.C")
                 if os.path.getmtime(os.path.join(R, "latest", f)) > os.path.getmtime(src)]
        with open(info, "w") as f:
            f.write(json.dumps({"commit": commit, "latest_modified": dirty, "sources_newer_than_binary": newer,
                                "built": time.strftime("%F %T", time.localtime(os.path.getmtime(src))),
                                "copied": time.strftime("%F %T")}) + "\n")
    return dst, json.loads(open(info).read())


# ---------------------------------------------------------------- nodes

def node_fields(text, n):
    """depth, inherited bound, predicted bound, fixes [(facility, location)] of a
    node line ('node <depth> <bound> <predicted> <nfix> <p[0..n-1]> <history>')"""
    t = text.split()
    p = [int(x) for x in t[5:5 + n]]
    return int(t[1]), float(t[2]), float(t[3]), [(i, j) for i, j in enumerate(p) if j >= 0]


def spec(fixes):
    return ",".join(f"{i}:{j}" for i, j in sorted(fixes))


def unit_hash(seed, fixes):
    """a uniform number in [0, 1) from the seed and the node's assignments: the
    stage-B sample is the frontier nodes with the smallest values (bottom-k), so
    it is fixed by the seed alone and any node's membership is known when the
    node is created"""
    h = hashlib.sha256(f"{seed}|{spec(fixes)}".encode()).digest()
    return int.from_bytes(h[:8], "big") / 2.0 ** 64


def strip_dnn(warm):
    """the same warm start without the parent's DNN state (the last field of a
    version-6 warm line becomes '-'); X, dnn_ext and the cut pool are kept"""
    if warm is None:
        return None
    head, _, last = warm.text.rpartition(" ")
    return warm if last == "-" else Warm(head + " -")


def write_one(src_path, index, out_path, warm="keep"):
    """a state file holding node `index` of src_path alone, statistics zeroed:
    warm = keep (its warm start as written), nodnn (without the DNN state) or
    none (no warm start at all, the node as -F starts it)"""
    st = State(src_path)
    # the warm line carries the DNN state from 6; 7 adds tables, 8 the settings record, 9 lets
    # the DNN state be a reference into a parent-state store ("@warm/<hash>.dnn") instead of
    # inline, which the solver writes only under QAP_DNN_WARM_DIR (these scripts never set it)
    if st.magic[1] not in ("6", "7", "8", "9"):
        raise ValueError(f"{src_path}: state version {st.magic[1]}, expected 6-9")
    st.zero_stats()
    text, w = st.nodes[index]
    if warm == "keep" and w is not None and w.text.rpartition(" ")[2].startswith("@"):
        raise ValueError(f"{src_path}: node {index}'s DNN state is a store reference; run without QAP_DNN_WARM_DIR")
    if warm == "nodnn":
        w = strip_dnn(w)
    elif warm == "none":
        w = None
    st.write(out_path, [(text, w)])
    return node_fields(text, st.n)


def slim_copy(src_path, out_path):
    """src_path with every warm start's DNN state removed (the frontier file)"""
    st = State(src_path)
    cache = {}
    nodes = []
    for text, w in st.nodes:
        if w is not None and id(w) not in cache:
            cache[id(w)] = strip_dnn(w)
        nodes.append((text, cache[id(w)] if w is not None else None))
    st.write(out_path, nodes)
    return st


# ---------------------------------------------------------------- solver runs

def run_bnbdrv(binary, args, env, log_path, trace_path=None, cpu=None, procs=None, tag=None):
    """bnbdrv ARGS in latest/ (where PARAM/ and SYM/ resolve), pinned to `cpu`,
    one thread (the settings set the thread variables), its own session; the DNN
    trace to trace_path. procs: a dict the caller uses to signal running solvers."""
    # pinned where taskset exists (Linux); macOS has no CPU pinning for user processes
    cmd = (["taskset", "-c", str(cpu)] if cpu is not None and shutil.which("taskset") else []) + [binary] + list(args)
    e = dict(env)
    if trace_path:
        e["QAP_DNN_TRACE"] = trace_path
        if os.path.exists(trace_path):
            os.remove(trace_path)
    t0 = time.time()
    with open(log_path, "w") as lf:
        lf.write("# " + " ".join(cmd) + "\n")
        lf.flush()
        p = subprocess.Popen(cmd, cwd=os.path.join(R, "latest"), env=e, stdout=lf, stderr=subprocess.STDOUT,
                             start_new_session=True)
        if procs is not None:
            procs[tag] = p
        try:
            rc = p.wait()
        finally:
            if procs is not None:
                procs.pop(tag, None)
    if trace_path and not os.path.exists(trace_path):
        open(trace_path, "w").close()           # no DNN bound ran (QPB fathomed, or depth 0)
    return rc, time.time() - t0


def log_checks(log_path):
    """what the solver log says about the run's setup"""
    txt = open(log_path, errors="replace").read()
    m = re.search(r"read symmetry info from (\S+)", txt)
    return {"sym_loaded": bool(m), "sym_file": m.group(1) if m else None,
            "sym_missing": "couldn't find" in txt,
            "dnn_line": next((l for l in txt.splitlines() if l.startswith("DNN bound at depths")), None),
            "resumed": next((l for l in txt.splitlines() if l.startswith("resumed:")), None),
            "new_incumbent": "New objective value" in txt or "new incumbent" in txt.lower()}


TRACE_RE = re.compile(
    r"pid (\d+) depth (\d+) fixes(?P<fixes>( \d+:\d+)*) \| first lb (?P<lb>\S+) primal (?P<primal>\S+) "
    r"status (?P<status>\S+) iters (?P<iters>\d+) cpu (?P<cpu>\S+) warm (?P<warm>\d) pool (?P<pool>\d) "
    r"kept (?P<kept>\d+) dropped (?P<dropped>\d+)(?P<pcp>[^|]*)\| cut (?P<cut>\d)(?P<cutrest>[^|]*)\| "
    r"cutoff (?P<cutoff>\S+) fathomed (?P<fathomed>\d) pool_out (?P<pool_out>\d+)")
CUT_RE = re.compile(r"best (?P<best>\S+) rounds (?P<rounds>\d+) iters (?P<iters>\d+) cpu (?P<cpu>\S+) "
                    r"status (?P<status>\S+)(?: cp (?P<cp>\d+) active (?P<active>\d+) stalled (?P<stalled>\d+))?")


def parse_trace(path):
    """one dict per DNN-bounded node in a QAP_DNN_TRACE file"""
    rows = []
    if not os.path.exists(path):
        return rows
    for line in open(path):
        m = TRACE_RE.search(line)
        if not m:
            continue
        g = m.groupdict()
        row = {"depth": int(m.group(2)), "fixes": spec([tuple(map(int, x.split(":"))) for x in g["fixes"].split()]),
               "lb": float(g["lb"]), "primal": float(g["primal"]), "status": g["status"], "iters": int(g["iters"]),
               "cpu": float(g["cpu"]), "warm": int(g["warm"]), "pool": int(g["pool"]), "kept": int(g["kept"]),
               "cut": int(g["cut"]), "cutoff": float(g["cutoff"]), "fathomed": int(g["fathomed"]),
               "pool_out": int(g["pool_out"])}
        if row["cut"]:
            c = CUT_RE.search(g["cutrest"])
            if c:
                row.update(cut_best=float(c["best"]), cut_rounds=int(c["rounds"]), cut_iters=int(c["iters"]),
                           cut_cpu=float(c["cpu"]), cut_status=c["status"],
                           cut_cp=int(c["cp"]) if c["cp"] else None)
        row["best"] = max(row["lb"], row.get("cut_best", -1e300))
        row["short"] = row["cutoff"] - row["best"]
        rows.append(row)
    return rows


def stats_of(st):
    """the per-task statistics kept in done.json (counters, per-level arrays, DNN)"""
    a = st.arrays
    return {"counters": st.counters, "times": st.times, "n_level": a["n_level"], "f_level": a["f_level"],
            "l_level": a["l_level"], "e_level": a["e_level"], "t_level": a["t_level"],
            "dnn_counts": a.get("dnn_counts", [0, 0]), "dnn_time": a.get("dnn_time", [0.0]),
            "obj": st.obj, "found": st.found}


def write_json(path, doc):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(doc, f, indent=1)
    os.replace(tmp, path)


def disk_free_gb(path):
    s = os.statvfs(path)
    return s.f_bavail * s.f_frsize / 1e9
