"""Read a DNN trace (QAP_DNN_TRACE, or the coordinator's trace.txt): one dict per DNN-bounded node."""
import re
R1 = re.compile(r'depth (\d+) fixes (.*?) \| first lb (\S+) primal \S+ status (\S+) iters (\d+) cpu (\S+) warm (\d) pool (\d) kept (\d+)')
R2 = re.compile(r'\| cut 1 best (\S+) rounds \d+ iters (\d+) cpu (\S+) status (\S+)')
R3 = re.compile(r'\| cutoff (\S+) fathomed (\d)')
def load(path):
    out = []
    for l in open(path):
        m = R1.search(l)
        if not m: continue
        e = dict(depth=int(m.group(1)), fixes=frozenset(m.group(2).split()), key=" ".join(sorted(m.group(2).split())),
                 lb=float(m.group(3)), status=m.group(4), iters=int(m.group(5)), cpu=float(m.group(6)),
                 warm=int(m.group(7)), pool=int(m.group(8)), kept=int(m.group(9)), cut=0, cit=0, ccpu=0.0, cstatus="", cbest=None, hist=[])
        c = R3.search(l); e["cutoff"] = float(c.group(1)); e["fathomed"] = int(c.group(2))
        p = R2.search(l)
        if p: e.update(cut=1, cbest=float(p.group(1)), cit=int(p.group(2)), ccpu=float(p.group(3)), cstatus=p.group(4))
        if "| hist" in l:
            e["hist"] = [(int(a), float(b)) for a, b in (t.split(":") for t in l.split("| hist")[1].split())]
        e["short"] = e["cutoff"] - e["lb"]
        e["left"] = (e["cutoff"] - e["cbest"]) if e["cut"] else e["short"]
        out.append(e)
    return out
def below(nodes):
    """for each node index: (descendant DNN nodes, their iterations) within this trace"""
    by = {}
    for i, e in enumerate(nodes): by[e["fixes"]] = i
    res = [[0, 0] for _ in nodes]
    for e in nodes:
        f = e["fixes"]
        # walk ancestors: subsets one smaller, repeatedly
        cur = [f]
        seen = set()
        while cur:
            nxt = []
            for s in cur:
                for x in s:
                    t = s - {x}
                    if t in seen: continue
                    seen.add(t)
                    if t in by:
                        res[by[t]][0] += 1; res[by[t]][1] += e["iters"] + e["cit"]
                        nxt.append(t)
            cur = nxt
    return res
