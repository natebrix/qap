#!/usr/bin/env python3
"""[m, from the production trace] What one level of branching is worth in bound units:
child DNN bound - parent DNN bound over the children of branched nodes of subtree 26:10,8:0.
A branched child's value is its give-up bound (exact); a fathomed child's is where its solve
stopped, at or just above the cutoff (a lower bound on its gain: censored)."""
import gzip, os, collections, statistics as st
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../.."))
CUT = 240514.05
rows = [l.split() for l in gzip.open(os.path.join(ROOT, "infra/phase1/runs/tho40ab-32/nodes.txt.gz"), "rt") if l.startswith("node ")]
by = {r[2]: r for r in rows}
G = collections.defaultdict(list)
for r in rows:
    par = by.get(r[2].rsplit(",", 1)[0])
    if par is None or par[7] != "B" or r[4] == "-": continue
    G[int(par[1])].append((float(r[4]) - float(par[4]), r[7] == "B", CUT - float(par[4])))
print("| parent depth | children | of them branched (uncensored) | median gain, all children (censored ones at their stop) | median gain of the branched children | quartiles of the branched children's gain | median parent shortfall |")
print("|---|---|---|---|---|---|---|")
for d in sorted(G):
    g = G[d]; unc = sorted(x[0] for x in g if x[1])
    q = lambda v, p: v[min(len(v) - 1, int(p * (len(v) - 1) + 0.5))]
    print("| %d | %d | %d | %.0f | %s | %s | %.0f |" % (d, len(g), len(unc), st.median(x[0] for x in g),
          "%.0f" % st.median(unc) if unc else "-", "%.0f / %.0f" % (q(unc, 0.25), q(unc, 0.75)) if unc else "-", st.median(x[2] for x in g)))
# parents short by more than 2,000: most children branch, so the gain is nearly uncensored
print()
for d in sorted(G):
    g = [x for x in G[d] if x[2] > 2000]
    if len(g) < 5: continue
    print("parents at depth %d short by > 2,000: %d children, %d%% fathomed, median gain %.0f (fathomed ones censored)" % (d, len(g), 100 * sum(not x[1] for x in g) / len(g), st.median(x[0] for x in g)))
