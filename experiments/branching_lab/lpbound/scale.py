"""scale.py NAME FIXES EPS...: build the AJ LP of one node, report its size, build time and RSS; then PDLP at each
EPS (warm-started from the previous), with the valid Lagrangian bound of its dual. Appends to out/scale.tsv."""
import sys, time, resource, numpy as np, ajlp
name, fx = sys.argv[1], sys.argv[2]; eps = [float(e) for e in sys.argv[3:]]
A, B, opt, perm = ajlp.inst(name)
Ar, Br, C, s, I, J = ajlp.reduce(A, B, fx)
t = time.time(); lp = ajlp.LP(Ar, Br, C); tb = time.time() - t
rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
print(f'{name} m={lp.m} vars={lp.nx+lp.ny} rows={lp.nrow} nnz={lp.M.nnz} build {tb:.1f}s rss {rss():.2f} GB', flush=True)
warm = None
with open('out/scale.tsv', 'a') as o:
    for e in eps:
        r = ajlp.pdlp_solve(lp, eps=e, warm=warm, time_limit=float(sys.argv[0] and 7200))
        warm = (r['x'], r['u'])
        line = f'{name}\t{fx}\t{lp.m}\tpdlp\t{e:g}\t{r["obj"]+s:.2f}\t{r["lagr"]+s:.2f}\t{r["iters"]}\t{r["wall"]:.1f}\t{r["cpu"]:.1f}\t{rss():.2f}\t{r["status"]}'
        print(line, flush=True); o.write(line + '\n'); o.flush()
