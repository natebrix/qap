"""List complete subtrees of the production trace by depth, size and root shortfall."""
import sys, collections
sys.path.insert(0, "/home/user/qap/experiments/tho40b_analysis")
import ntrace
t = ntrace.load("/home/user/qap/infra/phase1/runs/tho40ab-32/nodes.txt.gz")
n = len(t)
print("nodes", n, "cpu", sum(t.cpu))
# per depth summary
by = collections.defaultdict(list)
for i in range(n):
    by[t.depth[i]].append(i)
print("depth nodes branched  cpu  med_short_branched")
for d in sorted(by):
    idx = by[d]; br = [i for i in idx if t.outcome[i] == "B"]
    sh = sorted(t.cutoff[i] - t.dnn[i] for i in br)
    print(d, len(idx), len(br), round(sum(t.cpu[i] for i in idx)), sh[len(sh)//2] if sh else "-", sh[-1] if sh else "-")
lo, hi = int(sys.argv[1]), int(sys.argv[2])
print("subtrees with size in [%d,%d], depth 3-5:" % (lo, hi))
rows = []
for i in range(n):
    if t.outcome[i] == "B" and 3 <= t.depth[i] <= 5 and lo <= t.size[i] <= hi and t.complete[i]:
        rows.append((t.cutoff[i] - t.dnn[i], t.depth[i], t.size[i], t.sub_cpu[i], t.max_depth[i], ",".join(t.fixes[i])))
for r in sorted(rows):
    print("short %7.1f depth %d size %5d cpu %7.0f maxdepth %2d  %s" % r)
