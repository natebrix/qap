# Pool a tai50b estimator run (one directory per seed) from the 64-bit per-level
# tables in est.out, printed after every dive: dive k's share is k times the
# k-dive mean less k-1 times the previous one. Writes dives.csv, levels.csv and
# dive_levels.csv into the run directory and prints the summary.
# Usage: python3 -I pool.py RUN_DIR [SEED ...]
import re, math, statistics as st, csv, sys, os
run = sys.argv[1]; seeds = sys.argv[2:] or ['41', '42']
def blocks(path):
    t = open(path).read()
    parts = re.split(r"^.*\[ +\d+ dives .*$", t, flags=re.M)[1:]
    out = []
    for p in parts:
        m = re.search(r"^ +Nodes: (\d+) /", p, re.M)
        if not m: continue
        tm = float(re.search(r"^ +Total Time: ([\d.e+]+) seconds", p, re.M).group(1))
        lv = {}; f = False
        for line in p.splitlines():
            if line.startswith('level'): f = True; continue
            if f:
                q = line.split()
                if len(q) >= 5 and q[0].isdigit(): lv[int(q[0])] = (float(q[1]), float(q[4]))
                else: break
        out.append((float(m.group(1)), tm, lv))
    return out
rows, levrows = [], []
for s in seeds:
    bl = blocks(f'{run}/s{s}/est.out')
    depths = [int(l.split()[0]) for l in open(f'{run}/s{s}/dives.txt')]
    prev = (0.0, 0.0, {})
    for k, (n, t, lv) in enumerate(bl, 1):
        L = set(lv) | set(prev[2]); g = lambda d, l, i: d.get(l, (0, 0))[i]
        rows.append(dict(seed=s, dive=k, depth=depths[k-1] if k <= len(depths) else '', nodes=k*n-(k-1)*prev[0], cpu_s=k*t-(k-1)*prev[1]))
        for l in sorted(L):
            levrows.append(dict(seed=s, dive=k, level=l, nodes=k*g(lv,l,0)-(k-1)*g(prev[2],l,0),
                                cpu_s=k*g(lv,l,1)-(k-1)*g(prev[2],l,1)))
        prev = (n, t, lv)
N = len(rows); nn = [r['nodes'] for r in rows]; tt = [r['cpu_s'] for r in rows]
print(f"{run}: {N} dives; nodes {st.mean(nn):.3g} ± {st.stdev(nn)/math.sqrt(N):.2g}; CPU-s {st.mean(tt):.3g} ± {st.stdev(tt)/math.sqrt(N):.2g}"
      f" ({st.mean(tt)/1000/86400:.0f} days on 1,000 cores); median dive {st.median(nn):.3g} nodes, {st.median(tt):.3g} CPU-s;"
      f" top dive {100*max(tt)/sum(tt):.0f}% of CPU; CPU-s per node {sum(tt)/sum(nn):.1f}")
idx = {(r['seed'], r['dive']): i for i, r in enumerate(rows)}; lv_out = []
for l in sorted({r['level'] for r in levrows}):
    a = [0.0]*N; b = [0.0]*N
    for r in levrows:
        if r['level'] == l: i = idx[(r['seed'], r['dive'])]; a[i] = r['nodes']; b[i] = r['cpu_s']
    lv_out.append(dict(level=l, nodes_mean=st.mean(a), nodes_se=st.stdev(a)/math.sqrt(N), cpu_mean=st.mean(b),
                       cpu_se=st.stdev(b)/math.sqrt(N), dives_reaching=sum(1 for x in a if x > 0.5)))
print("level  nodes ± SE   CPU-s ± SE   reach")
for r in lv_out:
    print(f"{r['level']:>3}  {r['nodes_mean']:.3g} ± {r['nodes_se']:.2g}   {r['cpu_mean']:.3g} ± {r['cpu_se']:.2g}   {r['dives_reaching']}")
for name, data in (('dives.csv', rows), ('levels.csv', lv_out), ('dive_levels.csv', levrows)):
    with open(os.path.join(run, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
