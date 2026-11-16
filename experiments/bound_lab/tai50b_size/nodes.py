"""Per-node table from QAP_DNN_TRACE files:  python3 nodes.py [-md] [-hist] dnn.txt ...
Columns: depth, fixes, on the path of QAPSOLN/tai50b.sln?, first solve (iterations, bound, shortfall to the cutoff),
phase (iterations, best, how it ended), outcome, CPU-s. A shortfall < 0 is a fathom."""
import sys, re
BK = [5,23,14,12,21,4,29,40,13,8,16,19,25,18,46,11,15,37,17,30,47,20,33,26,36,48,22,10,3,31,2,7,0,32,6,42,24,28,44,27,39,35,43,45,1,34,49,38,9,41]
def parse(ln):
    m = re.match(r"pid (\d+) depth (\d+) fixes(.*?) \| first lb (\S+) primal (\S+) status (\S+) iters (\d+) cpu (\S+) warm (\d) pool (\d) kept (\d+) dropped (\d+)(?: pcp \d+ handed \d+)? \| cut (\d)(.*?) \| cutoff (\S+) fathomed (\d) pool_out (\d+)(?: \| hist (.*))?", ln)
    if not m: return None
    d = dict(pid=int(m[1]), depth=int(m[2]), fixes=m[3].split(), lb=float(m[4]), fstatus=m[6], fit=int(m[7]), fcpu=float(m[8]),
             warm=int(m[9]), pool=int(m[10]), kept=int(m[11]), cut=int(m[13]), cutoff=float(m[15]), fath=int(m[16]), hist=m[18])
    d["onpath"] = all(BK[int(x.split(":")[0])] == int(x.split(":")[1]) for x in d["fixes"])
    d["cit"] = d["ccpu"] = 0; d["best"] = None; d["cstatus"] = "-"
    if d["cut"]:
        c = re.match(r" best (\S+) rounds (\d+) iters (\d+) cpu (\S+) status (\S+)(?: cp (\d) active (\d+) stalled (\d) gaveup (\d))?", m[14])
        d["best"] = float(c[1]); d["cit"] = int(c[3]); d["ccpu"] = float(c[4]); d["cstatus"] = c[5]
        if c[6]: d["active"] = int(c[7]); d["cstatus"] += " stalled" * int(c[8]) + " gaveup" * int(c[9])
    d["final"] = max(d["lb"], d["best"] if d["best"] is not None else -1e300)
    d["out"] = ("D" if d["fstatus"] == "fathomed" else "C") if d["fath"] else "B"
    return d
def load(files):
    out = []
    for f in files:
        for ln in open(f):
            d = parse(ln)
            if d: d["file"] = f; out.append(d)
    return out
if __name__ == "__main__":
    a = sys.argv[1:]; hist = "-hist" in a; a = [x for x in a if x not in ("-md", "-hist")]
    print("| depth | fixes | on path | first solve: iterations | bound | short | phase iterations | phase best | short | phase end | outcome | CPU-s |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for d in load(a):
        print("| %d | `%s` | %s | %d | %s | %s | %s | %s | %s | %s | %s | %.0f |" % (d["depth"], ",".join(d["fixes"]) or "-", "yes" if d["onpath"] else "", d["fit"],
              f"{d['lb']:,.0f}", f"{d['cutoff']-d['lb']:,.0f}", d["cit"] or "", f"{d['best']:,.0f}" if d["cut"] else "",
              f"{d['cutoff']-d['best']:,.0f}" if d["cut"] else "", d["cstatus"], d["out"], d["fcpu"] + d["ccpu"]))
        if hist and d["hist"]:
            h = [x.split(":") for x in d["hist"].split()]
            print("  hist (short, M):", " ".join("%s:%.2f" % (i, (d["cutoff"] - float(v)) / 1e6) for i, v in h[:: max(1, len(h) // 30)]))
