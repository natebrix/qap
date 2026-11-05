#!/usr/bin/env python3
"""Phase-0 summaries (README.md).

  summarize.py instance OUT    on the instance: OUT/summary.json and a short text
  summarize.py table RUNDIR    on the laptop: one row per family from
                               RUNDIR/<family>/out/summary.json and RUNDIR/prices.tsv,
                               written to RUNDIR/summary.md

The reference for the fathom set is the 32-node run on the cloud box's 2.1 GHz
Xeon with the same flags (experiments/wil50_feasibility/speed/r32_noncert/
runs/base_<k>.log): 18 fathomed, 4,230 iterations to fathom. Iteration counts
can differ at the edges on another CPU (other BLAS kernels); the differences
are listed, not hidden.
"""
import glob
import json
import os
import re
import statistics
import sys

D = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(D, "..", "..", "experiments", "wil50_feasibility", "speed", "r32_noncert", "runs")


def result(path):
    try:
        for line in open(path, errors="replace"):
            if line.startswith("RESULT"):
                return {k: v for k, v in re.findall(r"(\w+)=(\S+)", line)}
    except OSError:
        pass
    return None


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split(" ", 1)
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def runs(d):
    rows = {}
    for p in glob.glob(os.path.join(d, "node_*.log")):
        k = int(re.search(r"node_(\d+)\.log$", p).group(1))
        r = result(p)
        rows[k] = r
    return rows


def instance(out):
    sysinfo = kv(os.path.join(out, "sysinfo.txt"))
    gates = {}
    if os.path.exists(os.path.join(out, "gates.txt")):
        for line in open(os.path.join(out, "gates.txt")):
            p = line.split()
            if len(p) >= 2:
                gates[p[0]] = p[1]
    s = dict(sysinfo=sysinfo, gates=gates)
    for mode in ("alone", "loaded"):
        rows = runs(os.path.join(out, mode))
        ok = {k: r for k, r in rows.items() if r}
        it = sum(int(r["iters"]) for r in ok.values())
        cpu = sum(float(r["cpu"]) for r in ok.values())
        s[mode] = dict(nodes=len(rows), finished=len(ok), iters=it, cpu_s=cpu,
                       ms_per_iter=1000 * cpu / it if it else None,
                       cpu_s_per_node=cpu / len(ok) if ok else None,
                       per_node={k: dict(status=r["status"], iters=int(r["iters"]), cpu=float(r["cpu"]),
                                         ms_per_iter=float(r["ms_per_iter"]), m=int(r["m"])) for k, r in ok.items()})
    # the fathom set against the reference
    diffs, fath, ref_fath = [], 0, 0
    for k, r in s["loaded"]["per_node"].items():
        ref = result(os.path.join(REF, f"base_{k}.log"))
        f = r["status"] == "fathomed"
        fath += f
        if ref:
            rf = ref["status"] == "fathomed"
            ref_fath += rf
            if rf != f or (f and int(ref["iters"]) != r["iters"]):
                diffs.append(f"node {k}: {ref['status']} {ref['iters']} -> {r['status']} {r['iters']}")
    s["fathom"] = dict(fathomed=fath, reference_fathomed=ref_fath, differences=diffs)
    a, l = s["alone"]["ms_per_iter"], s["loaded"]["ms_per_iter"]
    s["loaded_over_alone"] = l / a if a and l else None
    json.dump(s, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(f"{sysinfo.get('instance_type')} {sysinfo.get('model')} cores {sysinfo.get('physical_cores')} "
          f"amx {sysinfo.get('amx_bf16')}")
    print("gates " + " ".join(f"{g}:{r}" for g, r in gates.items()))
    for mode in ("alone", "loaded"):
        m = s[mode]
        if m["ms_per_iter"]:
            print(f"{mode}: {m['finished']}/{m['nodes']} nodes, {m['iters']} iterations, "
                  f"{m['ms_per_iter']:.1f} ms per iteration, {m['cpu_s_per_node']:.1f} CPU s per node")
    if s["loaded_over_alone"]:
        print(f"loaded / alone per iteration: {s['loaded_over_alone']:.3f}")
    print(f"fathomed {fath} (reference {ref_fath} on the same nodes); differences: {diffs or 'none'}")


def table(rundir):
    prices = {}
    pt = os.path.join(rundir, "prices.tsv")
    if os.path.exists(pt):
        for line in open(pt):
            p = line.split()
            if len(p) >= 3 and p[0] != "type":
                prices[p[0]] = (p[1], p[2])          # type -> (market, $/h at launch)
    rows = []
    for sj in sorted(glob.glob(os.path.join(rundir, "*", "out", "summary.json"))):
        s = json.load(open(sj))
        si = s["sysinfo"]
        typ = si.get("instance_type", "?")
        cores = int(si.get("physical_cores", 0) or 0)
        mk, price = prices.get(typ, ("?", "?"))
        try:
            ph = float(price)
        except ValueError:
            ph = None
        L = s["loaded"]
        core_h = ph / cores if ph and cores else None
        per_node = L["cpu_s_per_node"]
        usd_1k = per_node * 1000 * core_h / 3600 if per_node and core_h else None
        rows.append(dict(typ=typ, model=si.get("model", "?"), arch=si.get("arch", "?"), cores=cores,
                         amx=si.get("amx_bf16", "?"), alone=s["alone"]["ms_per_iter"], loaded=L["ms_per_iter"],
                         ratio=s.get("loaded_over_alone"), per_node=per_node, market=mk, ph=ph, core_h=core_h,
                         usd_1k=usd_1k, gates=s["gates"], fath=s["fathom"]))
    # steps of families that did not reach a summary (a failed build, an interruption)
    partial = []
    for st in sorted(glob.glob(os.path.join(rundir, "*", "out", "steps.txt"))):
        fam = st.split(os.sep)[-3]
        if not os.path.exists(os.path.join(os.path.dirname(st), "summary.json")):
            last = open(st).read().strip().splitlines()[-2:]
            partial.append(f"- {fam}: no summary; last steps: " + " | ".join(last))
    ref = min((r["loaded"] for r in rows if r["loaded"]), default=None)
    out = [f"# Phase 0, {os.path.basename(os.path.normpath(rundir))}", "",
           "| instance | CPU | arch | cores | AMX | ms/iteration alone | loaded | loaded / alone | CPU s per node (loaded) "
           "| market, $/h | $/core-hour | $ per 1,000 DNN nodes | gates | fathomed (reference) |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    f = lambda x, fmt: (fmt % x) if x is not None else "-"
    for r in rows:
        g = " ".join(f"{k}:{v}" for k, v in r["gates"].items() if v != "pass") or "all pass"
        fa = f"{r['fath']['fathomed']} ({r['fath']['reference_fathomed']})"
        if r["fath"]["differences"]:
            fa += f", {len(r['fath']['differences'])} differ"
        out.append(f"| {r['typ']} | {r['model']} | {r['arch']} | {r['cores']} | {r['amx']} | {f(r['alone'], '%.1f')} | "
                   f"{f(r['loaded'], '%.1f')} | {f(r['ratio'], '%.2f')} | {f(r['per_node'], '%.1f')} | "
                   f"{r['market']} {f(r['ph'], '%.4f')} | {f(r['core_h'], '%.4f')} | {f(r['usd_1k'], '%.3f')} | {g} | {fa} |")
    out.append("")
    if ref:
        out.append("Per-iteration speed under load, relative to the fastest: " +
                   ", ".join(f"{r['typ']} {r['loaded'] / ref:.2f}x" for r in rows if r["loaded"]))
    for r in rows:
        for d in r["fath"]["differences"]:
            out.append(f"- {r['typ']}: {d}")
    out += partial
    txt = "\n".join(out) + "\n"
    open(os.path.join(rundir, "summary.md"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("instance", "table"):
        sys.exit(__doc__)
    (instance if sys.argv[1] == "instance" else table)(sys.argv[2])
