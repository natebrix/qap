"""Where the production subtree's CPU sits, by depth and give-up shortfall of the branched nodes."""
import sys, collections
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
import ntrace
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
bins = [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3500, 5000, 1e9]
tot = sum(t.cpu)
print("branched nodes by depth and shortfall bin: count / share of all CPU in their children (direct children only) / share under them (whole subtree)")
for d in range(2, 8):
    idx = [i for i in range(len(t)) if t.depth[i] == d and t.outcome[i] == "B"]
    print("depth %d: %d branched, own CPU %.1f%% of all" % (d, len(idx), 100 * sum(t.cpu[i] for i in idx) / tot))
    for lo, hi in zip(bins, bins[1:]):
        b = [i for i in idx if lo <= t.cutoff[i] - t.dnn[i] < hi]
        if not b: continue
        kids = sum(t.cpu[k] for i in b for k in t.children[i])
        sub = sum(t.sub_cpu[i] - t.cpu[i] for i in b)
        allf = sum(1 for i in b if all(t.outcome[k] != "B" for k in t.children[i]))
        print("   %5.0f-%-5.0f n %4d (all children fathom: %3d)  children %5.1f%%  subtree below %5.1f%%  mean subtree nodes %7.1f" % (lo, min(hi, 99999), len(b), allf, 100 * kids / tot, 100 * sub / tot, sum(t.size[i] - 1 for i in b) / len(b)))
