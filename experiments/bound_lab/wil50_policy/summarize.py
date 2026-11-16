"""One row per subtree x arm from runs/NAME_ARM/{dnn.txt,log.txt}.
  python3 summarize.py [NAME ...]          the table (ratios against NAME_base)
  python3 summarize.py --phases NAME_ARM   every phase of one run
  python3 summarize.py --matched NAME ARM1 ARM2   node by node: the nodes both arms bounded, and the rest"""
import os, sys, glob, collections
from tracelib import load, below
H = os.path.dirname(os.path.abspath(__file__))
def stats(d):
    n = load(os.path.join(d, "dnn.txt")) if os.path.exists(os.path.join(d, "dnn.txt")) else []
    s = dict(done=os.path.exists(os.path.join(d, "done")), nodes=len(n))
    s["first"] = sum(e["iters"] for e in n); s["first_fath"] = sum(e["iters"] for e in n if e["status"] == "fathomed")
    s["n_first_fath"] = sum(1 for e in n if e["status"] == "fathomed")
    ph = [e for e in n if e["cut"]]
    s["phases"] = len(ph)
    for k in ("fathomed", "stalled", "iter_limit"):
        s["n_" + k] = sum(1 for e in ph if e["cstatus"] == k); s["it_" + k] = sum(e["cit"] for e in ph if e["cstatus"] == k)
    s["n_other"] = len(ph) - s["n_fathomed"] - s["n_stalled"] - s["n_iter_limit"]
    s["phase"] = sum(e["cit"] for e in ph); s["total"] = s["first"] + s["phase"]
    s["cpu"] = sum(e["cpu"] + e["ccpu"] for e in n)
    s["unfath"] = sum(1 for e in n if not e["fathomed"])
    s["nophase"] = sum(1 for e in n if e["status"] != "fathomed" and not e["cut"])
    s["depths"] = collections.Counter(e["depth"] for e in n)
    log = os.path.join(d, "log.txt"); t = open(log).read() if os.path.exists(log) else ""
    s["newobj"] = t.count("New objective value"); s["complete"] = "Fathoms by cause" in t
    cmd = os.path.join(d, "cmd.txt"); s["exit"] = [l.strip() for l in open(cmd) if l.startswith("exit")] if os.path.exists(cmd) else []
    return s, n
def table(names):
    runs = sorted(glob.glob(os.path.join(H, "runs", "*_*")))
    by = collections.defaultdict(dict)
    for d in runs:
        name, arm = os.path.basename(d).rsplit("_", 1)
        if names and name not in names: continue
        by[name][arm] = stats(d)[0]
    print("| subtree | arm | DNN nodes (by depth) | first solves: it (give-ups) | phases: n = fathomed + handed over/stalled + limit | phase it: fathomed + stalled + limit | total it | x base | CPU s | x base | new obj | complete |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name in by:
        b = by[name].get("base")
        for arm in sorted(by[name], key=lambda a: (a != "base", a)):
            s = by[name][arm]
            r = lambda k: ("%.2fx" % (b[k] / s[k])) if b and b["done"] and s["done"] and s[k] else "-"
            print("| %s | %s | %d (%s) | %s (%d) | %d = %d + %d + %d | %s + %s + %s = %s | %s | %s | %.0f | %s | %d | %s |" % (
                name, arm, s["nodes"], ", ".join("%d: %d" % kv for kv in sorted(s["depths"].items())), format(s["first"], ","), s["nodes"] - s["n_first_fath"],
                s["phases"], s["n_fathomed"], s["n_stalled"], s["n_iter_limit"], format(s["it_fathomed"], ","), format(s["it_stalled"], ","), format(s["it_iter_limit"], ","), format(s["phase"], ","),
                format(s["total"], ","), r("total"), s["cpu"], r("cpu"), s["newobj"], ("yes" if s["complete"] and s["done"] else "RUNNING" if not s["done"] else "NO")))
def phases(run):
    s, n = stats(os.path.join(H, "runs", run)); b = below(n)
    print("depth pool | first: it short | phase: it status left | below: nodes it")
    for e, (k, it) in zip(n, b):
        if e["status"] == "fathomed": continue
        print("%d %d | %4d %7.1f | %5d %-10s %7.1f | %4d %6d | %s" % (e["depth"], e["pool"], e["iters"], e["short"], e["cit"], e["cstatus"] or "-", e["left"], k, it, e["key"]))
def matched(name, a1, a2):
    n1 = {e["key"]: e for e in stats(os.path.join(H, "runs", name + "_" + a1))[1]}
    n2 = {e["key"]: e for e in stats(os.path.join(H, "runs", name + "_" + a2))[1]}
    both = [k for k in n1 if k in n2]
    tot = lambda n, ks: sum(n[k]["iters"] + n[k]["cit"] for k in ks)
    print("%s: %s %d nodes, %s %d nodes, both %d" % (name, a1, len(n1), a2, len(n2), len(both)))
    print("  matched nodes: %s %d it, %s %d it (%.2fx)" % (a1, tot(n1, both), a2, tot(n2, both), tot(n1, both) / max(1, tot(n2, both))))
    for cls, f in [("fathomed in the first solve in both", lambda k: n1[k]["status"] == "fathomed" and n2[k]["status"] == "fathomed"),
                   ("phase in both", lambda k: n1[k]["cut"] and n2[k]["cut"]),
                   ("other", lambda k: not (n1[k]["status"] == "fathomed" and n2[k]["status"] == "fathomed") and not (n1[k]["cut"] and n2[k]["cut"]))]:
        ks = [k for k in both if f(k)]
        if ks: print("    %-38s %4d nodes: %6d it against %6d (%.2fx)" % (cls, len(ks), tot(n1, ks), tot(n2, ks), tot(n1, ks) / max(1, tot(n2, ks))))
    o1 = [k for k in n1 if k not in n2]; o2 = [k for k in n2 if k not in n1]
    print("  only %s: %d nodes %d it; only %s: %d nodes %d it" % (a1, len(o1), tot(n1, o1), a2, len(o2), tot(n2, o2)))
if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--phases": phases(a[1])
    elif a and a[0] == "--matched": matched(a[1], a[2], a[3])
    else: table(a)
