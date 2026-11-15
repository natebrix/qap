"""The path check's verdict: in runs/pathK_ARM (the cold path's node at depth K and its children, incumbent 240,518)
the path node and its path child must both be branched (B), and no bound on either may pass 240,516.
  python3 pathverify.py [ARM]"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, H
arm = sys.argv[1] if len(sys.argv) > 1 else "z030"
rows = [l.rstrip("\n").split("\t") for l in open("/home/user/qap/experiments/tho40b_analysis/path_check/cold.tsv")][1:]
path = [r[1] for r in rows[1:]]          # the fixes in order
ok_all = True
for d in sorted(glob.glob(os.path.join(H, "runs", "path*_" + arm))):
    k = int(os.path.basename(d).split("_")[0][4:])
    nodes, dnn = load(d)
    want = {",".join(path[:k]): "node", ",".join(path[:k + 1]): "child"}
    done = os.path.exists(os.path.join(d, "done"))
    seen = {}
    for x, e in zip(nodes, dnn):
        if x["fixes"] in want:
            best = max(e["lb"], e["cbest"] if e["cbest"] is not None else -1e99)
            ok = x["out"] == "B" and best <= 240516.0 + 1e-6
            seen[want[x["fixes"]]] = ok
            ok_all &= ok
            print("depth %2d path %-5s %s first solve %.4f (%d it%s)%s  -> %s" % (
                x["depth"], want[x["fixes"]], x["out"], e["lb"], e["iters"], ", pooled %d cuts" % e["kept"] if e["kept"] else "",
                ", phase best %.4f after %d it (%s)" % (e["cbest"], e["cit"], e["cstatus"]) if e["cut"] else "", "ok" if ok else "FAIL"))
    others = [x for x in nodes if x["fixes"] not in want]
    if others or not done: print("   siblings of the path child: %d bounded, %d fathomed; run %s" % (len(others), sum(1 for x in others if x["out"] in "DC"), "complete" if done else "IN PROGRESS"))
    cmd = open(os.path.join(d, "cmd.txt")).read().split("\n")[0]
    node_only = (" -d %d " % k) in cmd      # run to the path node's own depth: its children were not bounded
    if node_only: print("   (the path node alone: -d %d)" % k)
    if done and len(seen) < (1 if node_only else 2): ok_all = False; print("   path node or path child missing from the trace: FAIL")
print("PASS" if ok_all else "FAIL")
