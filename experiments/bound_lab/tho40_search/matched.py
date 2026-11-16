"""Two runs of the same root, child by child: iterations under each of the root's children (the child and
everything below it), summed over the children whose subtrees are complete in both.
  python3 matched.py RUN_A RUN_B     e.g. w18_z030 w18_zall"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, H
def under(run):
    n, d = load(os.path.join(H, "runs", run))
    root = n[0]["fixes"]; k = root.count(",") + 2
    tot = {}; comp = {}
    for x, e in zip(n[1:], d[1:]):
        c = ",".join(x["fixes"].split(",")[:k])
        t = tot.setdefault(c, [0, 0, 0, 0]); t[0] += e["iters"] + e["cit"]; t[1] += 1; t[2] += e["cut"]
        if x["fixes"] == c: comp[c] = (x["out"], int(x["kids"]) if x["kids"] != "-" else 0)
        elif x["fixes"].count(",") == c.count(",") + 1: t[3] += 1
    done = {c for c in comp if comp[c][0] != "B" or tot[c][3] == comp[c][1]}   # (grandchildren that branch again: not checked)
    return n[0], d[0], tot, done
a, b = sys.argv[1], sys.argv[2]
ra, da, ta, ca = under(a); rb, db, tb, cb = under(b)
both = sorted(ca & cb)
print("root: %s %d + %d it; %s %d + %d it" % (a, da["iters"], da["cit"], b, db["iters"], db["cit"]))
print("children complete in both: %d (of %d in %s, %d in %s)" % (len(both), len(ca), a, len(cb), b))
sa = sum(ta[c][0] for c in both); sb = sum(tb[c][0] for c in both)
print("iterations under them: %s %d (%d nodes, %d phases), %s %d (%d nodes, %d phases): %s/%s = %.2f" % (
    a, sa, sum(ta[c][1] for c in both), sum(ta[c][2] for c in both), b, sb, sum(tb[c][1] for c in both), sum(tb[c][2] for c in both), b, a, sb / max(sa, 1)))
