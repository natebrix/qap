"""Summarize the runs: one row per subtree x arm from runs/NAME_ARM/{nodes.txt,dnn.txt,log.txt}.

  python3 summarize.py            table of every finished run (and partial counts of running ones)
  python3 summarize.py --depth NAME_ARM   per-depth detail of one run
"""
import os, re, sys, glob, collections
H = os.path.dirname(os.path.abspath(__file__))

def load(d):
    nodes = []
    p = os.path.join(d, "nodes.txt")
    if os.path.exists(p):
        for ln in open(p):
            f = ln.split()
            if len(f) >= 12 and f[0] == "node":
                nodes.append(dict(depth=int(f[1]), fixes=f[2], dnn=f[4], cut=f[5], cutoff=float(f[6]), out=f[7], kids=f[10], cpu=float(f[11])))
    dnn = []
    p = os.path.join(d, "dnn.txt")
    if os.path.exists(p):
        for ln in open(p):
            m = re.search(r"depth (\d+) fixes (.*?) \| first lb (\S+) primal \S+ status (\w+) iters (\d+) cpu (\S+) warm (\d) pool (\d) kept (\d+)", ln)
            if not m: continue
            e = dict(depth=int(m.group(1)), lb=float(m.group(3)), status=m.group(4), iters=int(m.group(5)), cpu=float(m.group(6)),
                     warm=int(m.group(7)), pool=int(m.group(8)), kept=int(m.group(9)), cut=0, cit=0, ccpu=0.0, cstatus="", cbest=None)
            c = re.search(r"\| cut 1 best (\S+) rounds (\d+) iters (\d+) cpu (\S+) status (\w+)", ln)
            if c:
                e.update(cut=1, cbest=float(c.group(1)), cit=int(c.group(3)), ccpu=float(c.group(4)), cstatus=c.group(5))
            e["cutoff"] = float(re.search(r"cutoff (\S+)", ln).group(1))
            e["fathomed"] = int(re.search(r"fathomed (\d)", ln).group(1))
            dnn.append(e)
    return nodes, dnn

def stats(d):
    nodes, dnn = load(d)
    s = collections.OrderedDict()
    s["done"] = os.path.exists(os.path.join(d, "done"))
    s["nodes"] = sum(1 for x in nodes if x["out"] != "P")
    s["B"] = sum(1 for x in nodes if x["out"] == "B")
    s["D"] = sum(1 for x in nodes if x["out"] == "D")
    s["C"] = sum(1 for x in nodes if x["out"] == "C")
    s["other"] = sum(1 for x in nodes if x["out"] not in "BDC")
    s["dnn_nodes"] = len(dnn)
    s["it_first"] = sum(e["iters"] for e in dnn)
    s["phases"] = sum(e["cut"] for e in dnn)
    s["phase_fath"] = sum(1 for e in dnn if e["cut"] and e["cstatus"] == "fathomed")
    s["it_cut"] = sum(e["cit"] for e in dnn)
    s["it_cut_fail"] = sum(e["cit"] for e in dnn if e["cut"] and e["cstatus"] != "fathomed")
    s["it_total"] = s["it_first"] + s["it_cut"]
    s["cpu"] = sum(x["cpu"] for x in nodes)
    s["cpu_dnn"] = sum(e["cpu"] for e in dnn); s["cpu_cut"] = sum(e["ccpu"] for e in dnn)
    root = dnn[0] if dnn else None
    s["root_it"] = (root["iters"] + root["cit"]) if root else 0
    s["root_first_it"] = root["iters"] if root else 0
    s["root_short"] = (root["cutoff"] - root["lb"]) if root else float("nan")
    s["root_cut"] = (root["cstatus"] or "-") if root else "-"
    s["root_cpu"] = nodes[0]["cpu"] if nodes and False else None
    log = os.path.join(d, "log.txt")
    s["newobj"] = 0
    if os.path.exists(log):
        t = open(log).read()
        s["newobj"] = len(re.findall(r"New objective value", t))
    s["maxdepth"] = max((x["depth"] for x in nodes), default=0)
    return s, nodes, dnn

def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--depth":
        s, nodes, dnn = stats(os.path.join(H, "runs", sys.argv[2]))
        by = collections.defaultdict(lambda: collections.Counter())
        for x in nodes:
            by[x["depth"]][x["out"]] += 1; by[x["depth"]]["cpu"] += x["cpu"]
        for e in dnn:
            b = by[e["depth"]]
            b["it_first"] += e["iters"]; b["phases"] += e["cut"]; b["it_cut"] += e["cit"]
            b["pf"] += 1 if e["cstatus"] == "fathomed" else 0
            b["pooled"] += 1 if e["kept"] > 0 else 0
        print("depth nodes  B    D    C  phases(fath) pooled it_first it_cut   cpu")
        for d in sorted(by):
            b = by[d]
            print("%4d %6d %4d %4d %4d %5d(%4d) %6d %8d %6d %7.0f" % (d, b["B"] + b["D"] + b["C"] + b["E"] + b["Q"], b["B"], b["D"], b["C"], b["phases"], b["pf"], b["pooled"], b["it_first"], b["it_cut"], b["cpu"]))
        print("phases: first-solve shortfall -> best, iters, status")
        for e in dnn:
            if e["cut"]:
                print("  d%d short %7.1f -> %8.1f  it %4d  %-9s first %d it%s" % (e["depth"], e["cutoff"] - e["lb"], e["cutoff"] - e["cbest"], e["cit"], e["cstatus"], e["iters"], " pooled" if e["kept"] else ""))
        return
    rows = []
    for d in sorted(glob.glob(os.path.join(H, "runs", "*_*"))):
        s, _, _ = stats(d)
        rows.append((os.path.basename(d), s))
    hdr = "run            done nodes    B    D    C  dnn_nodes it_first phases fath it_cut it_cut_fail it_total   cpu_s  root_it root_short root_cut newobj"
    print(hdr)
    for name, s in rows:
        print("%-14s %-4s %5d %4d %4d %4d %9d %8d %6d %4d %6d %11d %8d %7.0f %8d %10.1f %8s %d" % (
            name, "y" if s["done"] else "run", s["nodes"], s["B"], s["D"], s["C"], s["dnn_nodes"], s["it_first"], s["phases"], s["phase_fath"],
            s["it_cut"], s["it_cut_fail"], s["it_total"], s["cpu"], s["root_it"], s["root_short"], s["root_cut"], s["newobj"]))

if __name__ == "__main__":
    main()
