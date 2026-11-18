import sys, math, statistics
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import replay
T = replay.Tree(replay.load_trace(sys.argv[1]) if sys.argv[1].endswith(".gz") else replay.load_shards(sys.argv[1]))
def rank(xs):
    o = sorted(range(len(xs)), key=lambda i: xs[i]); r = [0]*len(xs)
    for k, i in enumerate(o): r[i] = k
    return r
def spear(a, b):
    ra, rb = rank(a), rank(b); n = len(a)
    if n < 3: return None
    ma, mb = statistics.fmean(ra), statistics.fmean(rb)
    num = sum((x-ma)*(y-mb) for x, y in zip(ra, rb)); den = math.sqrt(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))
    return num/den if den else None
top1 = tot = 0; sp = []; mass = []
for i in range(len(T.depth)):
    if not T.branched[i]: continue
    surv = [j for j in T.kids[i] if T.branched[j]]
    if len(surv) < 2: continue
    sh = [T.cutoff[j] - T.best[j] for j in surv]; sz = [T.C[j] for j in surv]
    tot += 1; top1 += sh.index(max(sh)) == sz.index(max(sz))
    s = spear(sh, sz)
    if s is not None: sp.append(s)
    mass.append(max(sz) / sum(sz))
print("nodes with >= 2 surviving children:", tot, "; the largest-shortfall survivor has the largest subtree:", f"{top1/tot:.0%}",
      "; Spearman(shortfall, subtree CPU) mean", round(statistics.fmean(sp), 2), "; the largest subtree's share of the survivors' CPU, median", round(statistics.median(mass), 2))
