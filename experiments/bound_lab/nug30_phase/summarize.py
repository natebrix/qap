"""Tables from runs/ARM/{nodes.txt,dnn.txt,log.txt,time.txt}.
  python3 summarize.py                 one row per arm (ratios against base)
  python3 summarize.py --depth ARM     per-depth detail and every phase of one arm
  python3 summarize.py --phases        phases by first-solve shortfall, all cut arms
"""
import os, re, sys, glob, collections
H = os.path.dirname(os.path.abspath(__file__))
ORDER = ["base", "base2", "z003", "z003b", "z002", "z010", "z030", "z060", "sko", "d1", "g1", "g2", "g2b", "g4", "g5", "g5w100", "g5w144", "g5w288", "g5e", "basea", "skob", "g2a", "g2s"]

def load(arm):
    d = os.path.join(H, "runs", arm)
    nodes, dnn = [], []
    for ln in open(os.path.join(d, "nodes.txt")):
        f = ln.split()
        if len(f) >= 12 and f[0] == "node":
            nodes.append(dict(depth=int(f[1]), fixes=f[2], dnn=f[4], cut=f[5], out=f[7], cpu=float(f[11])))
    for ln in open(os.path.join(d, "dnn.txt")):
        m = re.search(r"depth (\d+) fixes (.*?) \| first lb (\S+) primal \S+ status (\w+) iters (\d+) cpu (\S+) warm (\d) pool (\d) kept (\d+)", ln)
        if not m: continue
        e = dict(depth=int(m.group(1)), fixes=m.group(2), lb=float(m.group(3)), status=m.group(4), iters=int(m.group(5)), cpu=float(m.group(6)),
                 warm=int(m.group(7)), kept=int(m.group(9)), cut=0, cit=0, ccpu=0.0, cstatus="", cbest=None)
        c = re.search(r"\| cut 1 best (\S+) rounds (\d+) iters (\d+) cpu (\S+) status (\w+)(?: cp \d+ active (\d+) stalled (\d) gaveup (\d))?", ln)
        if c:
            e.update(cut=1, cbest=float(c.group(1)), cit=int(c.group(3)), ccpu=float(c.group(4)), cstatus=c.group(5), stalled=c.group(7))
        e["cutoff"] = float(re.search(r"cutoff (\S+)", ln).group(1))
        e["fathomed"] = int(re.search(r"fathomed (\d)", ln).group(1))
        dnn.append(e)
    return nodes, dnn

def stats(arm):
    d = os.path.join(H, "runs", arm)
    nodes, dnn = load(arm)
    s = {}
    for k in "BDCQE": s[k] = sum(1 for x in nodes if x["out"] == k)
    s["nodes"] = sum(1 for x in nodes if x["out"] != "P")
    s["dnn_nodes"] = len(dnn)
    s["it_first"] = sum(e["iters"] for e in dnn)
    s["phases"] = sum(e["cut"] for e in dnn)
    s["pf"] = sum(1 for e in dnn if e["cstatus"] == "fathomed")
    s["it_cut"] = sum(e["cit"] for e in dnn)
    s["it"] = s["it_first"] + s["it_cut"]
    s["cpu_first"] = sum(e["cpu"] for e in dnn); s["cpu_cut"] = sum(e["ccpu"] for e in dnn)
    s["pooled"] = sum(1 for e in dnn if e["kept"] > 0)
    t = open(os.path.join(d, "log.txt")).read()
    s["newobj"] = len(re.findall(r"New objective value", t))
    m = re.search(r"wall ([\d.]+) user ([\d.]+) sys ([\d.]+)", open(os.path.join(d, "time.txt")).read()) if os.path.exists(os.path.join(d, "time.txt")) else None
    s["wall"] = float(m.group(1)) if m else float("nan")
    s["user"] = float(m.group(2)) + float(m.group(3)) if m else float("nan")
    m = re.search(r"nodes\s*[:=]?\s*(\d+)", t.split("---- RESULTS ----")[-1]) if "RESULTS" in t else None
    s["log_tail"] = t.split("---- RESULTS ----")[-1] if "RESULTS" in t else ""
    s["cpu_trace"] = sum(x["cpu"] for x in nodes)
    s["done"] = os.path.exists(os.path.join(d, "done"))
    return s, nodes, dnn

def arms():
    a = [os.path.basename(p) for p in glob.glob(os.path.join(H, "runs", "*")) if os.path.exists(os.path.join(p, "dnn.txt"))]
    return [x for x in ORDER if x in a] + sorted(x for x in a if x not in ORDER)

def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--depth":
        s, nodes, dnn = stats(sys.argv[2])
        by = collections.defaultdict(collections.Counter)
        for x in nodes:
            by[x["depth"]][x["out"]] += 1; by[x["depth"]]["cpu"] += x["cpu"]
        for e in dnn:
            b = by[e["depth"]]
            b["it_first"] += e["iters"]; b["phases"] += e["cut"]; b["it_cut"] += e["cit"]
            b["pf"] += e["cstatus"] == "fathomed"; b["pooled"] += e["kept"] > 0
        print("| depth | nodes | branched | fathomed: QPB / first solve / phase | phases (fathomed) | pooled first solves | iterations: first + phases | CPU s |\n|---|---|---|---|---|---|---|---|")
        for d in sorted(by):
            b = by[d]
            print("| %d | %d | %d | %d / %d / %d | %d (%d) | %d | %s + %s | %.0f |" % (d, sum(b[k] for k in "BDCQE"), b["B"], b["Q"] + b["E"], b["D"], b["C"], b["phases"], b["pf"], b["pooled"], f'{b["it_first"]:,}', f'{b["it_cut"]:,}', b["cpu"]))
        return
    if len(sys.argv) > 1 and sys.argv[1] == "--phases":
        bins = [0, 5, 10, 20, 40, 61, 100, 184, 400]
        for arm in arms():
            s, nodes, dnn = stats(arm)
            if not s["phases"]: continue
            print(f"\n{arm}: shortfall bin | phases | fathomed | median it fathomed | median it missed | median left after a miss | pooled")
            for lo, hi in zip(bins, bins[1:]):
                ph = [e for e in dnn if e["cut"] and lo <= e["cutoff"] - e["lb"] < hi]
                if not ph: continue
                med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
                f = [e for e in ph if e["cstatus"] == "fathomed"]; mi = [e for e in ph if e["cstatus"] != "fathomed"]
                print("  %3d-%3d | %4d | %4d | %5.0f | %5.0f | %6.1f | %d" % (lo, hi, len(ph), len(f), med([e["cit"] for e in f]), med([e["cit"] for e in mi]), med([e["cutoff"] - e["cbest"] for e in mi]), sum(e["kept"] > 0 for e in ph)))
        return
    base = None
    print("| arm | nodes | branched | fathomed: QPB / first solve / phase | phases (fathomed) | iterations: first solves + phases = total | phase share | CPU s (trace) | wall s | iterations vs base | CPU vs base | nodes vs base | new obj |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for arm in arms():
        s, _, _ = stats(arm)
        if arm == "base": base = s
        r = lambda k: ("%.2fx" % (base[k] / s[k])) if base and arm != "base" else ""
        print("| %s%s | %s | %d | %d / %d / %d | %d (%d) | %s + %s = %s | %.0f%% | %.0f | %.0f | %s | %s | %s | %d |" % (
            arm, "" if s["done"] else " (running)", f'{s["nodes"]:,}', s["B"], s["Q"] + s["E"], s["D"], s["C"], s["phases"], s["pf"],
            f'{s["it_first"]:,}', f'{s["it_cut"]:,}', f'{s["it"]:,}', 100.0 * s["it_cut"] / max(1, s["it"]), s["cpu_trace"], s["wall"], r("it"), r("cpu_trace"), r("nodes"), s["newobj"]))

if __name__ == "__main__":
    main()
